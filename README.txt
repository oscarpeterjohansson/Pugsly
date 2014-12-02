A PubChem Tool (PCT): The Power User Gateway (PUG)


    "The Pug is a breed of dog with a wrinkly, 
    short-muzzled face and curled tail"


This module was written to provide a programmatic option to download an sdf
file for specified UIDs. Perhaps, it can be used for other purposes as well.

Choose a database (db) (class str)
      pccompound   -  PubChem Compound
      pcsubstance  -  PubChem Substance

Choose a format (form) (class str)
      text-asn     -  full records
                      textual ASN.1
      binary-asn   -  binary ASN.1
      xml          -  textual XML
      sdf          -  SD file format, for chemical structures
      image        -  images, format is always .zip containing 
                      multiple .png
                      full-size depiction
      image-small  -  thumbnail depiction
      smiles       -  selected string fields, format is: SID/CID <tab> <string>
                      Isomeric SMILES
      inchi        -  InChI

Choose a compression type (compr) (class str)
      none         -  no compression
                      (to use this option you might have to 
                      remove the "PCT-Download_compression"
                      element from the tup object in "constrxlm_init_input".
                      At least it did not suffice to set the value to "")
      gzip         -  gzip format
      bzip2        -  bzip2 format

Specify UIDs (User Identifiers) (python iterable, class int or str)
                                (OR
                                 comma-separated list (python str))
      eg. [2244, 5362129]     -  Aspirin, Ramipril
      OR  "2244,5362129"
          


For more info, see
https://pubchem.ncbi.nlm.nih.gov/pug/pughelp.html  

For info on PUG, alt. see
https://pubchem.ncbi.nlm.nih.gov/pc_fetch/pc_fetch.cgi  
for a non-programmatic interface


For use with the Python intrepreter,
   Simple usage example:

    pctpug.retr_data(
        uid=[2244,123908],
        form="text-asn", 
        compr="gzip", 
        db="pccompound",
        filename=""

For command-line use, 
  try:

    --help
