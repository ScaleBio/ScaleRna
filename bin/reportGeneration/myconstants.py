
#!/usr/bin/env python
# Plotly express figure style to be used for all figures 
DEFAULT_FIGURE_STYLE="none"
# Number of cells sampled for scatter plots 
SAMPLING_NUMBER = 4000
# Color mapping used for qc filter categorized scatter plots 
QC_COLORMAP = {True: 'rgb(39, 139, 176)', False: 'rgb(233, 237, 245)'}
# Color mapping used for barnyard species categorization scatter plots (human and mouse hardcoded as it is the only barnyard genome used [could be changed])
BARNYARD_COLORMAP = {"None": 'rgb(179, 188, 201)', "Ambig": 'rgb(201, 147, 166)', 'Mixed': 'rgb(242, 5, 33)', 'Human': 'rgb(36, 36, 227)', 'Mouse': 'rgb(27, 99, 25)'}
# 
UNKNOWN_COLORMAP={"Unknown": 'rgb(179, 188, 201)'}

BARCODE_SHORTHAND_TO_NAME={ 'drop':'Droplet Barcodes', 'P7':'P7 Barcodes', 'lig':'Ligation Barcodes','rt':'Reverse Transcription Barcodes', 'umi':'UMI'}

DEFAULT_COLOR_SEQUENCE=['#012B64','#013780', '#00419A', '#014BB0','#0454C4','#045BD3','#0663E3','#086AF2','#2D80F4','#448FF7','#5FA0F9','#81B3F8','#9EC4F9','#BBD5F9']

# Columns to display in stats tables 
DISPLAY_COLUMNS=["Metric","Value"]