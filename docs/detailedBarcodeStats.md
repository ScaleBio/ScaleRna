## Read Level metrics
**Pass**: Reads for which all expected barcodes were found \
**LinkerError**: Reads which were filtered because the fixed linker sequence between the barcodes could not be found \
**BarcodeError**: Reads for which at least one barcode could not be matched against the expected sequences (whitelist) \
**SequenceError**: Reads excluded from barcode matching, e.g. because they were too short

## Barcode level  metrics
**Exact**: The read contains an exact (no errors) match against one of the expected sequences (whitelist) for this barcode \
**Corrected**: The barcode sequence contains at least one mismatch, but could be corrected to an unique whitelist sequence \
**Ambiguous**: The barcode sequence has the same number of mismatches to two different whitelist sequences and is hence filtered \
**NoMatch**: The barcode sequence cannot be matched to any whitelist sequence \
**Error**: No barcode sequence can be extracted; typically because the linker sequence used to locate it in the read was not found \

The UMI is a random sequence with no whitelist. In that case all sequences containing 'N's or pure homopolymers are **Filtered**, other sequences are **Pass**.