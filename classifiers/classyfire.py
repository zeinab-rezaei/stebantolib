
import json
import time
from pathlib import Path
from typing import Dict, Optional
from urllib.parse import quote

import requests


# Constants
CLASSYFIRE_CACHE_PATH = "classyfire_cache.json"
CLASSYFIRE_ENTITY_URL = "https://classyfire.wishartlab.ca/entities/{inchikey}.json"
CLASSYFIRE_CLASSIFY_URL = "https://classyfire.wishartlab.ca/classification.json"
NPCLASSIFIER_URL = "https://npclassifier.ucsd.edu/classify"
FIEHNLAB_MIRROR = "https://cfb.fiehnlab.ucdavis.edu/entities/{inchikey}.json"


class ClassifierConfig:
    """Configuration for classification services"""

    def __init__(self, timeout: int = 12, debug: bool = False):
        self.timeout = timeout
        self.debug = debug
        self.disable_wishart = False


class ClassificationCache:
    """Handles caching of classification results"""

    def __init__(self, cache_path: str = CLASSYFIRE_CACHE_PATH):
        self.cache_path = Path(cache_path)
        self._cache = self._load_cache()

    def _load_cache(self) -> Dict:
        """Load cache from disk"""
        try:
            return json.loads(self.cache_path.read_text(encoding="utf-8"))
        except (FileNotFoundError, json.JSONDecodeError):
            return {}

    def save(self):
        """Persist cache to disk"""
        self.cache_path.write_text(
            json.dumps(self._cache, ensure_ascii=False, indent=2),
            encoding="utf-8"
        )

    def get(self, key: str) -> Optional[Dict]:
        return self._cache.get(key)

    def put(self, key: str, value: Dict):
        self._cache[key] = value
        self.save()


class ClassificationService:
    """Main service for compound classification"""

    def __init__(self, config: ClassifierConfig):
        self.config = config
        self.cache = ClassificationCache()
        self.session = self._create_session()

    def _create_session(self) -> requests.Session:
        """Create configured HTTP session"""
        session = requests.Session()
        session.headers.update({"User-Agent": "massbank-msp-converter/1.0"})
        return session

    def _log(self, msg: str):
        """Debug logging"""
        if self.config.debug:
            print(f"[classify] {msg}")

    def get_compound_classes(self, inchikey: Optional[str], smiles: Optional[str]) -> Optional[Dict]:
        """
        Get compound classification with fallback strategy:
        1. ClassyFire by InChIKey
        2. ClassyFire by SMILES
        3. NPClassifier by SMILES
        """
        inchikey = (inchikey or "").strip()
        smiles = (smiles or "").strip()

        # Try ClassyFire by InChIKey first
        if inchikey:
            result = self._try_classyfire_by_inchikey(inchikey)
            if result:
                return result

        # Fallback to SMILES-based classification
        if smiles:
            result = self._try_smiles_classification(smiles)
            if result:
                return result

        return None

    def _try_classyfire_by_inchikey(self, inchikey: str) -> Optional[Dict]:
        """Try to get classification from ClassyFire using InChIKey"""
        cache_key = f"ik:{inchikey}"
        cached = self.cache.get(cache_key)
        if self._is_valid_classification(cached):
            self._log(f"Cache hit (InChIKey): {inchikey}")
            return cached

        result = self._query_classyfire_entities(inchikey)
        if self._is_valid_classification(result):
            self.cache.put(cache_key, result)
            return result

        self._log(f"No ChemOnt via InChIKey {inchikey}")
        return None

    def _try_smiles_classification(self, smiles: str) -> Optional[Dict]:
        """Try SMILES-based classification with fallback providers"""
        cache_key = f"smiles:{smiles}"
        cached = self.cache.get(cache_key)
        if self._is_valid_classification(cached):
            self._log("Cache hit (SMILES)")
            return cached

        # Try ClassyFire by SMILES
        if not self.config.disable_wishart:
            result = self._query_classyfire_by_smiles(smiles)
            if self._is_valid_classification(result):
                self.cache.put(cache_key, result)
                return result

        # Fallback to NPClassifier
        return self._query_npclassifier(smiles)

    def _query_classyfire_entities(self, inchikey: str) -> Optional[Dict]:
        """Query ClassyFire entity endpoints with mirror fallback"""
        mirrors = []
        if not self.config.disable_wishart:
            mirrors.append(CLASSYFIRE_ENTITY_URL.format(inchikey=quote(inchikey)))
        mirrors.append(FIEHNLAB_MIRROR.format(inchikey=quote(inchikey)))

        for url in mirrors:
            try:
                self._log(f"GET {url}")
                response = self._request_with_retry(url, timeout=self.config.timeout)
                if response and response.ok:
                    result = self._normalize_chemont_entity(response.json())
                    if self._is_valid_classification(result):
                        return result
            except requests.exceptions.RequestException as e:
                self._handle_request_error(e, url)

        return None

    def _query_classyfire_by_smiles(self, smiles: str) -> Optional[Dict]:
        """Query ClassyFire classification endpoint with SMILES"""
        if self.config.disable_wishart:
            return None

        try:
            self._log("POST ClassyFire classifier")
            response = self._request_with_retry(
                CLASSYFIRE_CLASSIFY_URL,
                method="POST",
                json={"smiles": smiles},
                timeout=max(self.config.timeout, 20)
            )
            if response and response.ok:
                return self._normalize_chemont_entity(response.json())
        except requests.exceptions.RequestException as e:
            self._handle_request_error(e, CLASSYFIRE_CLASSIFY_URL)

        return None

    def _query_npclassifier(self, smiles: str) -> Optional[Dict]:
        """Query NPClassifier service"""
        url = f"{NPCLASSIFIER_URL}?smiles={quote(smiles)}"
        try:
            self._log("GET NPClassifier")
            response = self.session.get(url, timeout=max(self.config.timeout, 20))
            if response.ok:
                npc_data = response.json()
                return self._normalize_npclassifier(npc_data)
        except requests.exceptions.RequestException as e:
            self._log(f"NPClassifier error: {e}")

        return None

    def _request_with_retry(self, url: str, method: str = "GET", **kwargs) -> Optional[requests.Response]:
        """Make HTTP request with retry logic for rate limiting"""
        try:
            if method.upper() == "POST":
                response = self.session.post(url, **kwargs)
            else:
                response = self.session.get(url, **kwargs)

            if response.status_code in (429, 503):
                time.sleep(0.7)  # Brief pause before retry
                if method.upper() == "POST":
                    response = self.session.post(url, **kwargs)
                else:
                    response = self.session.get(url, **kwargs)

            return response
        except requests.exceptions.RequestException:
            return None

    def _handle_request_error(self, error: Exception, url: str):
        """Handle HTTP request errors, particularly DNS issues"""
        error_str = str(error).lower()
        if ("getaddrinfo failed" in error_str or "name resolution" in error_str) and "wishartlab" in url:
            self.config.disable_wishart = True
            self._log("Wishart host disabled for this run (DNS error)")
        self._log(f"Request error: {error}")

    @staticmethod
    def _is_valid_classification(result: Optional[Dict]) -> bool:
        """Check if classification result is valid"""
        if not result:
            return False
        required_fields = ("kingdom", "superclass", "class", "subclass", "direct_parent")
        return any(result.get(field) for field in required_fields)

    @staticmethod
    def _normalize_chemont_entity(entity_data: Dict) -> Optional[Dict]:
        """Normalize ClassyFire/ChemOnt entity data"""
        if not entity_data:
            return None

        def get_node(name: str) -> str:
            return (entity_data.get(name) or {}).get("name") or ""

        def get_node_id(name: str) -> str:
            return (entity_data.get(name) or {}).get("chemont_id") or ""

        result = {
            "ontology": "ChemOnt",
            "kingdom": get_node("kingdom"),
            "superclass": get_node("superclass"),
            "class": get_node("class"),
            "subclass": get_node("subclass"),
            "direct_parent": get_node("direct_parent"),
            "ontology_id": get_node_id("direct_parent") or get_node_id("class") or get_node_id("superclass"),
            "source": "ClassyFire",
        }

        return result if ClassificationService._is_valid_classification(result) else None

    @staticmethod
    def _normalize_npclassifier(npc_data: Dict) -> Optional[Dict]:
        """Normalize NPClassifier data"""
        if not npc_data:
            return None

        return {
            "ontology": "NPClassifier",
            "kingdom": "",
            "superclass": (npc_data.get("superclass") or [""])[0],
            "class": (npc_data.get("class") or [""])[0],
            "subclass": "",
            "direct_parent": (npc_data.get("pathway") or [""])[0],
            "ontology_id": "",
            "source": "NPClassifier",
        }


def get_compound_classes(inchikey: Optional[str], smiles: Optional[str], debug: bool = False) -> Optional[Dict]:
    """
    Compatibility function that matches your existing interface.
    This wraps the new ClassificationService for backward compatibility.
    """
    config = ClassifierConfig(debug=debug)
    service = ClassificationService(config)
    return service.get_compound_classes(inchikey, smiles)