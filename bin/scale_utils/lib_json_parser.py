import json
from collections import OrderedDict
from pathlib import Path

class LibJsonParser:
    """
    Parsed library structure json information
    """

    def __init__(self, fname: Path):
        self.parent_dir = fname.parent
        self.json_contents = self.readJSON(fname)
        self.all_whitelist_contents_by_alias, self.all_whitelist_contents_by_read = self.parse_all_whitelists()
        self.sample_barcode_fname = self.get_sample_barcode_fname()

    def get_sample_barcode_fname(self) -> Path:
        """
        Function to get the sample barcode file name from the JSON contents

        Returns:
            Path to the sample barcode file
        """
        for barcode in self.json_contents["barcodes"]:
            if barcode["name"] == self.json_contents["sample_barcode"]:
                return self.parent_dir / barcode["sequences"]
        raise ValueError("Sample barcode not found in JSON contents")
    
    def getMaxWellNumberAndLetter(self, fname: Path) -> tuple[str, int]:
        """
        Get maximum well coordinate from whitelist file

        Args:
            fname: Path to file

        Returns:
            Maximum letter and maximum number corresponding to the max well coordinate
        """
        max_letter = "A"
        max_number = 1
        with open(fname) as f:
            for line in f:
                line = line.strip()
                split_line = line.split("\t")
                letter = split_line[0][-1]
                numbers = int(split_line[0][:-1])
                max_letter = letter if letter > max_letter else max_letter
                max_number = numbers if numbers > max_number else max_number
        return max_letter, max_number

    def get_barcodes_range_of_whitelist_file(self, fname: Path, sep: str = "\t"):
        """
        Get the range of barcodes (first and last entry) in the whitelist file
        Assumption: The first column contains the alias

        Args:
            fname: Path to file
            sep: Separator used in the file
        
        Returns:
            Tuple with first and last entry in the whitelist file
        """
        with open(fname) as f:
            lines = f.readlines()
            first_entry = lines[0].strip().split(sep)[0]
            last_entry = lines[-1].strip().split(sep)[0]
        return first_entry, last_entry
    
    def load_whitelist(self, fname: Path):
        """
        Function to load a whitelist file and return a dictionary with alias and sequences

        Args:
            fname: Path to file

        Returns:
            Dictionary with alias and sequences
        """
        bcs = {}
        for line in open(fname):
            if line.startswith("#"):
                continue
            line = line.strip().split()
            # No name given, use sequence as name
            if len(line) == 1:
                seq = line[0]
                name = seq
            else:
                name = line[0]
                seq = line[1:]
            if name in bcs:
                raise ValueError(f"Duplicate barcode alias {name} found in whitelist file {fname}")
            if seq in bcs.values():
                raise ValueError(f"Duplicate barcode sequence {seq} found in whitelist file {fname}")
            bcs[name] = seq
        return bcs

    def parse_all_whitelists(self) -> dict[str, dict]:
        """
        Function to parse all whitelists from the JSON contents

        Args:
            json_contents: Dictionary with contents of json file

        Returns:
            Dictionary with all whitelists parsed from the JSON contents
        """
        all_whitelist_contents_by_alias = {}
        all_whitelist_contents_by_read = {}
        for bc in self.json_contents["barcodes"]:
            if bc.get("type", None) not in ["library_index", "umi", "target"]:
                    bcs = self.load_whitelist(self.parent_dir / bc["sequences"])
                    if len(bcs) != len(set(bcs)):
                        raise ValueError(f"Duplicate barcodes found in whitelist file: {bc['sequences']}")
                    all_whitelist_contents_by_alias[bc["alias"]] = bcs
                    all_whitelist_contents_by_read[bc["read"]] = bcs
        return all_whitelist_contents_by_alias, all_whitelist_contents_by_read

    def readJSON(self, file: Path):
        """
        Function to read in JSON file and return it as a dictionary

        Args:
            file: Path to .json file

        Returns:
            Dictionary with contents of json file
        """
        with open(file) as f:
            str = f.read()
            strStripped = str.rstrip()
            pairs_hook = OrderedDict
            parsedJSON = json.loads(strStripped, object_pairs_hook=pairs_hook)
        return parsedJSON


