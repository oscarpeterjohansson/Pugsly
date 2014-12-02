#!/usr/bin/python
"""
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

"""

__version__ = "2.0"
__email__ = "oscarpeterjohansson@outlook.com"
__author__ = "Johansson, O."
__contributors__ = ""
__licence__ = "GPL-3"

import argparse
import os
import pycurl
import re
import StringIO
import sys
import time
import traceback
import xml.etree.ElementTree as ET


CGI_URL = "http://pubchem.ncbi.nlm.nih.gov/pug/pug.cgi"           # Common Gateway Interface
DTD_URL = "https://pubchem.ncbi.nlm.nih.gov/pug/pug.xsd"          # Document Type Definition
XML_DECLARATION = """<?xml version="1.0" encoding="utf-8"?>\n"""
# !DOCTYPE PCT-Data defines that the root element is PCT-Data
DOCTYPE = """<!DOCTYPE PCT-Data PUBLIC "-//NCBI//NCBI PCTools/EN" "http://pubchem.ncbi.nlm.nih.gov/pug/pug.dtd">\n""" 
INIT_XML_STRUCT = (                               # don't touch
    ("PCT-Data_input",
     ("PCT-InputData",
      ("PCT-InputData_download",
       ("PCT-Download",
        (
            ("PCT-Download_uids",
             ("PCT-QueryUids",
              ("PCT-QueryUids_ids",
               ("PCT-ID-List",
                (
                    ("PCT-ID-List_db",),
                    ("PCT-ID-List_uids",)))))),
            ("PCT-Download_format",),
            ("PCT-Download_compression",)))))))
POLL_XML_STRUCT = (                               # don't touch
    ("PCT-Data_input",
     ("PCT-InputData",
      ("PCT-InputData_request",
       ("PCT-Request",
        (
            ("PCT-Request_reqid",),
            ("PCT-Request_type",)))))))


def indent(elem,level=0):
    """
    Prettyprint:  This (particular function) is not one of my inventions. 
    But, found it inspiring for my own recursive function (apndr).
    """
    i = "\n" + level * "  "
    if len(elem):
        if not elem.text or not elem.text.strip():
            elem.text = i + "  "
        if not elem.tail or not elem.tail.strip():
            elem.tail = i
        for elem in elem:
            indent(elem, level+1)                 # recursion
        if not elem.tail or not elem.tail.strip():
            elem.tail = i
    else:
        if level and (not elem.tail or not elem.tail.strip()):
            elem.tail = i


def apndr(parent, tup, i=0):
    """
    Append with recursion  -  My own invention!
    parent                 -  <class 'xml.etree.ElementTree.Element'>
    tup                    -  tuple  -  of length two (or more): See examples ...
                                        the level of "nestedness" 
                                        corresponds to the level in the xml tree
    """
    i += 1                                              #  to inform where in ... the mess, that we are
    #print "i:",i
    if isinstance(tup,tuple) and len(tup) == 1:         # the end/leaf
        if isinstance(tup[0],basestring):                               
            ET.SubElement(parent, tup[0])               # cannot do more on this path
    elif isinstance(tup,tuple) and len(tup) > 1:        # not leaf yet
        if isinstance(tup[0],basestring):               # we have a child to append
            ET.SubElement(parent, tup[0])               # append child ...
            n_parent = parent.find(".//%s" % tup[0])
            #apndr(parent[0], tup[1], i)                # DOES NOT WORK, because idx0 refers to ... pretty much anything
            apndr(n_parent, tup[1], i)                  # ... then move over to next tuple; update parent
        if isinstance(tup[0],tuple):
            for tup in tup:
                apndr(parent, tup, i)                   # same parent for all tuples
    else:                                               # just, a way to detect errors in the larger tuple
        print "isinstance(tup,tuple) == %s " % (isinstance(tup,tuple),)
        print "len(tup) == %s " % (len(tup),)
        raise StandardError("Error in structure of tuple tup?")
    return parent


def mk_tree(root_tag="PCT-Data",tup=INIT_XML_STRUCT):
    """
    Construct tree with recursion.  Text and attrib will have to be set
    afterwards.
    root_tag  -  str    -  tag of root
    tup       -  tuple  -  nested tuple structure 
    """
    root = ET.Element(root_tag)
    root = apndr(root, tup)
    indent(root)
    tree = ET.ElementTree(root)
    tree.write(sys.stdout)
    return tree


def http_post(url, xml):
    """
    Post request to PUG with  pycurl: The reason for trying pycurl instead
    of "urllib2" or "requests" is that "Name or service not known" error 
    was issued with those.
    """
    buf = StringIO.StringIO() 
    c = pycurl.Curl()
    c.setopt(pycurl.URL, url)
    c.setopt(pycurl.POST, 1)
    c.setopt(pycurl.HTTPHEADER, ["Content-type: text/xml"])
    c.setopt(pycurl.POSTFIELDS, xml)
    c.setopt(pycurl.WRITEFUNCTION, buf.write)
    c.setopt(pycurl.FOLLOWLOCATION, 1)
    c.setopt(pycurl.USERAGENT, "Mozilla/5.0")
    try:
        c.perform()
        buf.seek(0,0)
        return buf.read()
    except (StandardError, pycurl.error) as e:
        tb = sys.exc_info()[2]
        traceback.print_tb(tb, limit=5, file=sys.stdout)
        print e
        return None


def mk_init_input_xml(uid=[1,99], form="sdf", compr="gzip", db="pccompound"):
    """ Construct "init input" xml for PC (PubChem) PUG - All included

    compr  -  str   -  compression (eg. gzip)
    db     -  str   -  database (eg. pccompound)
    form   -  str   -  format of output (eg. sdf)
    uid    -  list  -  of str UIDs (User Identifiers)
    """
    root = ET.Element("PCT-Data")
    root = apndr(root, INIT_XML_STRUCT)         # construct
    parent = root.find(".//%s" % ("PCT-ID-List_uids",))
    for i in uid:
        elem = ET.SubElement(parent,"PCT-ID-List_uids_E")
        elem.text = str(i)                      # UIDs
    parent = root.find(".//%s" % ("PCT-ID-List_db",)) 
    parent.text = str(db)                       # database
    parent = root.find(".//%s" % ("PCT-Download_format",))      
    parent.set("value",str(form))               # format
    parent = root.find(".//%s" % ("PCT-Download_compression",)) 
    parent.set("value",str(compr))              # compression
    indent(root)                                # prettyprint
    tree = ET.ElementTree(root)
    f = StringIO.StringIO()                     # convert (instance of) tree to string
    tree.write(f)                                                 
    f.seek(0)
    tree_str = f.read()
    f.close
    tree_str = XML_DECLARATION + DOCTYPE + tree_str
    print tree_str                              # display (alt: tree.write(sys.stdout) )
    return tree_str


def mk_poll_xml(reqid=402936103567975582):
    """ Construct "poll" xml for PC (PubChem) PUG

    reqid  -  str  -  request id, obatined with waiting message
    """
    root = ET.Element("PCT-Data")
    root = apndr(root, POLL_XML_STRUCT)         # construct
    parent = root.find(".//%s" % ("PCT-Request_reqid",))        
    parent.text = str(reqid)                    # database
    parent = root.find(".//%s" % ("PCT-Request_type",))         
    parent.set("value", "status")               # format
    indent(root)                                # prettyprint
    tree = ET.ElementTree(root)
    fd = StringIO.StringIO()                    # convert (instance of) tree to string
    tree.write(fd)
    fd.seek(0,0)
    tree_str = fd.read()
    fd.close()
    tree_str = XML_DECLARATION + DOCTYPE + tree_str
    print tree_str                              # display (alt: tree.write(sys.stdout) )
    return tree_str
            

def str2lst(uid):
    """ comma-separated list (str) to python list (list) """
    tmp = uid.split(',')
    tmp = [u.strip() for u in tmp]
    return filter(lambda x: x != "", tmp)


def dl_request(uid=[1,99], form="sdf", compr="gzip", db="pccompound"):
    """
    Send download request. Download data from database for UIDs in format
    with compression. Retrieve dl URL
    
    compr  -  str   -  compression (eg. gzip)
    db     -  str   -  database (eg. pccompound)
    form   -  str   -  format of output (eg. sdf)
    uid    -  list  -  of str UIDs (User Identifiers)
                       OR
           -  str   -  a comma-separated list (python str)
    """
    if isinstance(uid,str): 
        uid = str2lst(uid)                     # convert comma-sep list to python list
    xii = mk_init_input_xml(uid, form, compr, db)
    rsp = http_post(CGI_URL, xii)              # (first) response
    #print rsp
    if rsp is None:
        return None
    rsp_intp = ET.fromstring(rsp)              # interpreted
    s = rsp_intp.find(".//%s" % ("PCT-Status",)) 
    if s is None:                              # status
        sys.stdout.write("\rStatus missing. That's strange ... why?\r\n")
        return None
    try:
        s_v = s.get("value")
        print "Status value: %s\n" % (s_v,)
        s_msgs = rsp_intp.findall(".//%s" % ("PCT-Status-Message_messages_E",))
        for i in s_msgs:
            print i.text
    except (AttributeError,StandardError) as e:
        tb = sys.exc_info()[2]
        traceback.print_tb(tb, limit=5, file=sys.stdout)
    poll = 0
    etim = 0
    wait = rsp_intp.find(".//%s" % ("PCT-Waiting_reqid",))
    while wait is not None:                    # poll
        print "Poll no.:", poll
        print "Elapsed time:", etim
        poll += 1
        reqid = wait.text
        xp = mk_poll_xml(reqid)
        rsp = http_post(CGI_URL, xp)
        rsp_intp = ET.fromstring(rsp)                            
        wait = rsp_intp.find(".//%s" % ("PCT-Waiting_reqid",))
        if wait is not None:
            delay = 5
            time.sleep(delay)
            etim += delay
    dlurl_elm = rsp_intp.find(".//%s" % ("PCT-Download-URL_url",))    
    if dlurl_elm is not None:                  # download url available?
        return dlurl_elm.text
    else:
        sys.stdout.write("\r\nExpected to find \"PCT-Download-URL_url\" in rsp, but couldn't\n")
        return None


def progress(download_t, download_d, upload_t, upload_d):
    """ Print progress of download/upload. For use with pycurl """
    stat = """
    \rTotal to download %d, Total downloaded %d, Total to upload %d, Total uploaded %d 
    """ % (download_t, download_d, upload_t, upload_d)
    stat = re.sub("[\n]","",stat) # "\n's are introduced with use of '"""' and newlines
    sys.stdout.write(stat)
    sys.stdout.flush()


def retr_data(uid=[1,99], form="sdf", compr="gzip", db="pccompound", output=""):
    """
    Retrieve data; write to file. If output is "", then os.getcwd and 
    basename(URL) will be used to create a name
    """
    dlurl = dl_request(uid, form, compr, db)
    print "Download-URL ",dlurl
    if not dlurl:
        return None
    if output is "":
        output = os.path.join( os.getcwd(), os.path.basename(dlurl) )
    fd = open(output,'w')
    c = pycurl.Curl()
    c.setopt(pycurl.URL, dlurl)
    c.setopt(pycurl.WRITEFUNCTION, fd.write)
    c.setopt(pycurl.NOPROGRESS, 0)
    c.setopt(pycurl.PROGRESSFUNCTION, progress)
    try:
        c.perform()
    except (StandardError, pycurl.error) as e:
        tb = sys.exc_info()[2]
        traceback.print_tb(tb, limit=5, file=sys.stdout)
        print e
    finally:
        fd.close()
        sys.stdout.write("\r\n")
        return None


parser = argparse.ArgumentParser(
    prog = sys.argv[0],               # .py-file
    description = """
    PubChem Power User Gateway (PUG): Fetch PubChem data using Unique 
    IDentifiers (UIDs)
    """,
    prefix_chars = "-",
    conflict_handler = "error",
    add_help = True
)

parser.add_argument(
    "--input",     
    dest = 'uid',
    nargs = 1,
    type = str,
    required = True,
    help = """
    Path to file containing comma-separated UIDs (unique identifiers, 
    ie. CIDs or SIDs). No default.
    """,
    metavar = 'UID'
)

parser.add_argument(
    "--form",
    dest = 'form',
    nargs = 1,
    default = "sdf",
    type = str,
    choices = ["text-asn","binary-asn","xml","sdf","image","image-small",
               "smiles","inchi"],
    required = True,
    help = """ 
    Choose a format:
    "text-asn"         (full records textual ASN.1),
    "binary-asn"       (binary ASN.1),
    "xml"              (textual XML),
    "sdf"              (SD file format, for chemical structures),
    "image"            (images, format is always ".zip" containing 
    multiple ".png", full-size depiction),
    "image-small"      (thumbnail depiction),
    "smiles"           (selected string fields, format is: SID/CID
                        <tab> <string> Isomeric SMILES),
    "inchi"            (InChI).
    Defaults to "sdf". 
    """,
    metavar = 'FORMAT'
)

parser.add_argument(
    "--compr",
    dest = 'compr',
    nargs = 1,
    default = "gzip",
    type = str,
    choices = ["none","gzip","bzip2"], 
    required = True,
    help = """
    Choose a compression type:
    "none"            (no compression),
    "gzip"            (gzip format),
    "bzip2"           (bzip2 format).
    Defaults to "gzip".
    """,
    metavar = 'COMPRESSION'
)

parser.add_argument(
    "--db",   
    dest = 'db',
    nargs = 1,
    default = "pccompound",
    type = str,
    choices = ["pccompound","pcsubstance"],
    required = True,
    help = """
    Choose a database, one of:
    "pccompound"      (PubChem Compound),
    "pcsubstance"     (PubChem Substance).
    Defaults to "pccompound".
    """,
    metavar = 'DATABASE'
)

parser.add_argument(
    "--output",
    dest = 'output',
    nargs = 1,
    default = "",
    type = str,
    required = True,
    help = """
    Path to file for storage of fetched data. Defaults to current working
    directory.
    """ ,
    metavar = 'OUTPUT'
)


def main(argv):
    """ For command-line use
    """
    ## Parse command-line arguments
    argv = parser.parse_args(args = argv)
    if os.path.isdir(argv.output[0]):
        ms = "--output OUTPUT: expected path to file; got path to directory"
        parser.error(ms) # sys.exit: status code 2
    uid = None
    try:    # a try to be able to use "parser.error(ms)"
        fp = open(argv.uid[0],'r')
        tmp = fp.read()
        uid = re.sub("\s", "", tmp) # "\s"=="[ \t\n\r\f\v]"
        fp.close()
    except (IOError,) as e:
        ms = "%s:\r\n%s\r\nerrno: %s\r\n" % (e.strerror, e.filename, e.errno)
        parser.error(ms)
    except (TypeError,) as e:
        ms = "%s" % (e.ms,)
        parser.error(ms) # sys.exit: status code 2
    output = None
    if isinstance(argv.output,list):
        output = argv.output[0]
    else:
        output = argv.output
    retr_data(
        uid = uid, 
        form = argv.form[0], 
        compr = argv.compr[0],
        db = argv.db[0],
        output = output
        )
    parser.exit(status=0, message=None)
    

if __name__ == "__main__":
    main(sys.argv[1:])
    
