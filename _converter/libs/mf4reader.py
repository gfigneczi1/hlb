#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec  4 19:02:24 2020

@author: ote2bp
"""

import os
import enum
import struct
import logging

TestMode = False

#Data reader class by type
class DataType(enum.IntEnum):
    # ASAMMDF documentation "4.2 Definition of Data Types Used And Mapping to ASAM Data Types"
    UINT8 = 1
    BYTE = 1
    CHAR = 1
    INT8 = 1
    UINT16 = 2
    INT16 = 2
    UINT32 = 4
    INT32 = 4
    UINT64 = 8
    INT64 = 8
    FLOAT32 = 4
    FLOAT64 = 8
    REAL = 8
    LINK = 8

    def readIntLittle(data_bytes, PositionInfo):
        return int.from_bytes(data_bytes[PositionInfo[0] : (PositionInfo[0] + PositionInfo[1])], byteorder='little', signed=True)

    def readFloat64(data_bytes, PositionInfo):
        return struct.unpack('d', data_bytes[PositionInfo[0] : (PositionInfo[0] + PositionInfo[1])])[0]

    def readStr(data_bytes, PositionInfo):
        return data_bytes[PositionInfo[0] : (PositionInfo[0] + PositionInfo[1])].decode("UTF8").strip()

    def readChecksum(data_bytes, PositionInfo):
        return data_bytes[PositionInfo[0] : (PositionInfo[0] + PositionInfo[1])]

class mf4reader:
    length_IDBLOCK = 64
    length_HeaderSection = 24
    length_ToLengthInHeader = 8
    length_BYTE = 1
    length_UINT64 = 8
    length_LINK = 8

    positions_Field = {}
    info_Field = {
        'Link': {'LINK': ['INT64', 1]},

        ## ASAMMDF documentation "4.3.1 Header section"
        'Header': {
            'id'        : ['CHAR',   4],
            'reserved'  : ['BYTE',   4],
            'length'    : ['UINT64', 1],
            'link_count': ['UINT64', 1],
        },

        ## ASAMMDF documentation "4.5.1 Block Structure of IDBLOCK"
        'IDBLOCK': {
            'id_file'              : ['CHAR',   8],
            'id_vers'              : ['CHAR',   8],
            'id_prog'              : ['CHAR',   8],
            'id_reserved1'         : ['BYTE',   4],
            'id_ver'               : ['UINT16', 1],
            'id_reserved2'         : ['BYTE',  30],
            'id_unfin_flags'       : ['UINT16', 1],
            'id_custom_unfin_flags': ['UINT16', 1],
        },

        ## ASAMMDF documentation "4.6.1 Block Structure of HDBLOCK"
        'HDBLOCK': {
            'hd_start_time_ns'   : ['UINT64', 1],
            'hd_tz_offset_min'   : ['INT16',  1],
            'hd_dst_offset_min'  : ['INT16',  1],
            'hd_time_flags'      : ['UINT8',  1],
            'hd_time_class'      : ['UINT8',  1],
            'hd_flags'           : ['UINT8',  1],
            'hd_reserved'        : ['BYTE',   1],
            'hd_start_angle_rad' : ['REAL',   1],
            'hd_start_distance_m': ['REAL',   1],
        },

        ## ASAMMDF documentation "4.7.1 Block Structure of MDBLOCK"
        # 'MDBLOCK': {
        #     'md_data': ['BYTE', 0],
        # },

        ## ASAMMDF documentation "4.8.1 Block Structure of TXBLOCK"
        # 'TXBLOCK': {
        #     'tx_data': ['BYTE', 0],
        # },

        ## ASAMMDF documentation "4.9.1 Block Structure of FHBLOCK"
        'FHBLOCK': {
            'fh_time_ns'       : ['UINT64', 1],
            'fh_tz_offset_min' : ['INT16',  1],
            'fh_dst_offset_min': ['INT16',  1],
            'fh_time_flags'    : ['UINT8',  1],
            'fh_reserved'      : ['BYTE',   3],
        },

        ## ASAMMDF documentation "4.10.1 Block Structure of CHBLOCK"
        'CHBLOCK': {
            'ch_element_count' : ['UINT32', 1],
            'ch_type'          : ['UINT8' , 1],
            'ch_reserved'      : ['BYTE'  , 3],
        },

        ## ASAMMDF documentation "4.11.1 Block Structure of ATBLOCK"
        'ATBLOCK': {
            'at_flags':         ['UINT16', 1],
            'at_creator_index': ['UINT16', 1],
            'at_reserved':      ['BYTE',   4],
            'at_md5_checksum':  ['BYTE',   16],
            'at_original_size': ['UINT64', 1],
            'at_embedded_size': ['UINT64', 1],
            # SPECIAL 'at_embedded_data': ['BYTE',   0],
        },

        ## ASAMMDF documentation "4.12.1 Block structure of EVBLOCK"
        'EVBLOCK': {
            'ev_type':             ['UINT8',  1],
            'ev_sync_type':        ['UINT8',  1],
            'ev_range_type':       ['UINT8',  1],
            'ev_cause':            ['UINT8',  1],
            'ev_flags':            ['UINT8',  1],
            'ev_reserved':         ['BYTE',   3],
            'ev_scope_count':      ['UINT32', 1],
            'ev_attachment_count': ['UINT16', 1],
            'ev_creator_index':    ['UINT16', 1],
            'ev_sync_base_value':  ['INT64',  1],
            'ev_sync_factor':      ['REAL',   1],
        },

        ## ASAMMDF documentation "4.13.1 Block Structure of DGBLOCK"
        'DGBLOCK': {
            'dg_rec_id_size': ['UINT8', 1],
            'dg_reserved':    ['BYTE',  7],
        },

        ## ASAMMDF documentation "4.14.1 Block Structure of CGBLOCK"
        'CGBLOCK': {
            'cg_record_id':      ['UINT64', 1],
            'cg_cycle_count':    ['UINT64', 1],
            'cg_flags':          ['UINT16', 1],
            'cg_path_separator': ['UINT16', 1],
            'cg_reserved':       ['BYTE',   4],
            'cg_data_bytes':     ['UINT32', 1],
            'cg_inval_bytes':    ['UINT32', 1],
        },

        ## ASAMMDF documentation "4.15.1 Block Structure of SIBLOCK"
        'SIBLOCK': {
            'si_type':     ['UINT8', 1],
            'si_bus_type': ['UINT8', 1],
            'si_flags':    ['UINT8', 1],
            'si_reserved': ['BYTE',  5],
        },

        ## ASAMMDF documentation "4.16.1 Block Structure of CNBLOCK"
        'CNBLOCK': {
            'cn_type':             ['UINT8',  1],
            'cn_sync_type':        ['UINT8',  1],
            'cn_data_type':        ['UINT8',  1],
            'cn_bit_offset':       ['UINT8',  1],
            'cn_byte_offset':      ['UINT32', 1],
            'cn_bit_count':        ['UINT32', 1],
            'cn_flags':            ['UINT32', 1],
            'cn_inval_bit_pos':    ['UINT32', 1],
            'cn_precision':        ['UINT8',  1],
            'cn_reserved':         ['BYTE',   1],
            'cn_attachment_count': ['UINT16', 1],
            'cn_val_range_min':    ['REAL',   1],
            'cn_val_range_max':    ['REAL',   1],
            'cn_limit_min':        ['REAL',   1],
            'cn_limit_max':        ['REAL',   1],
            'cn_limit_ext_min':    ['REAL',   1],
            'cn_limit_ext_max':    ['REAL',   1],
        },

        ## ASAMMDF documentation "4.17.1 Block Structure of CCBLOCK"
        'CCBLOCK': {
            'cc_type':          ['UINT8',  1],
            'cc_precision':     ['UINT8',  1],
            'cc_flags':         ['UINT16', 1],
            'cc_ref_count':     ['UINT16', 1],
            'cc_val_count':     ['UINT16', 1],
            'cc_phy_range_min': ['REAL',   1],
            'cc_phy_range_max': ['REAL',   1],
            'cc_val':           ['REAL',   0],
        },
    }


    """
        Constructor, opens the mf4, checks it, initializes inner structures
    """
    def __init__(self, filepath=None):
        #first Meta after ##HD
        self.header = {}
        #All ##FH blocks
        self.comment = ''
        self.filehistory = {}
        self.hierarchies = {}
        self.attachments = {}
        self.events = {}
        self.datagroups = {}

        self.links_Channels = []
        self.links_Channelgroups = []
        self.links_Datagroups = []

        self.defineFieldPositions()

        self.ChannelsList = []
        self.ChannelsInfo = []
        self.link_DGBLOCK_current = 0
        self.link_CGBLOCK_current = 0
        self.link_CNBLOCK_current = 0


        # Opens and maps the file
        if filepath is not None:
            self.open(filepath)


    # Opens a file, called by the constructor
    def open(self, filepath):
        self.filepath = filepath
        self.filesize = os.stat(filepath).st_size
        #TODO: larger buffering
        self.filepointer = open(self.filepath, "rb", buffering=4096)
        # Reads the ID block from beginning of the file
        self.fileinfo = self.blockID()
        if self.fileinfo is not None:
            # Reads the header block, it contains info about the measurement
            self.header = self.blockHeaDer()
            if self.header is not None:
                """
                    lets create the map of the file following the links from the header block
                    It results 5 link chains
                """
                links_HDBLOCK = self.header['links']

                # Comment
                if links_HDBLOCK['hd_md_comment'] != 0:
                    comment_raw = self.getComment(links_HDBLOCK['hd_md_comment'])
                    lines_comment = []
                    while 0 != comment_raw.count('\n'):
                        index_LineEnd = comment_raw.index('\n')
                        line_comment = comment_raw[:index_LineEnd+1]
                        lines_comment.append(line_comment)
                        comment_raw = comment_raw[index_LineEnd+1:]
                    lines_comment.append(comment_raw)
                    self.comment = lines_comment

                # First history block
                if links_HDBLOCK['hd_fh_first'] != 0:
                    self.filehistory = self.readLinks(links_HDBLOCK['hd_fh_first'], "##FH")
                else:
                    logging.error("ERROR3 Missing history block")

                # First data group block
                if links_HDBLOCK['hd_dg_first'] != 0:
                    self.mapDataGroups(links_HDBLOCK['hd_dg_first'])
                else:
                    logging.error('The first datagroup block is missing from the header!')

                # First channel hierarchy block
                if links_HDBLOCK['hd_ch_first'] != 0:
                    self.hierarchies = self.processChannelHierarchy(links_HDBLOCK['hd_ch_first'])

                # First attachment block
                if links_HDBLOCK['hd_at_first'] != 0:
                    self.attachments = self.readLinks(links_HDBLOCK['hd_at_first'], "##AT")

                # First event block
                if links_HDBLOCK['hd_ev_first'] != 0:
                    self.events = self.readLinks(links_HDBLOCK['hd_ev_first'], "##EV")

        return

    def getFieldLength(self, Block, Field):
        FieldType = self.info_Field[Block][Field][0]
        FieldCount = self.info_Field[Block][Field][1]
        if 0 == FieldCount:
            FieldCount = 1

        if 'UINT8' == FieldType:
            TypeLength = 1
        elif 'BYTE' == FieldType:
            TypeLength = 1
        elif 'CHAR' == FieldType:
            TypeLength = 1
        elif 'INT8' == FieldType:
            TypeLength = 1
        elif 'UINT16' == FieldType:
            TypeLength = 2
        elif 'INT16' == FieldType:
            TypeLength = 2
        elif 'UINT32' == FieldType:
            TypeLength = 4
        elif 'INT32' == FieldType:
            TypeLength = 4
        elif 'UINT64' == FieldType:
            TypeLength = 8
        elif 'INT64' == FieldType:
            TypeLength = 8
        elif 'FLOAT32' == FieldType:
            TypeLength = 4
        elif 'FLOAT64' == FieldType:
            TypeLength = 8
        elif 'REAL' == FieldType:
            TypeLength = 8
        elif 'LINK' == FieldType:
            TypeLength = 8

        FieldLength = TypeLength * FieldCount
        return FieldLength

    def getUnpackFormat(self, Type):
        if 'UINT8' == Type:
            Format = '<B'
        elif 'INT8' == Type:
            Format = '<b'
        elif 'UINT16' == Type:
            Format = '<H'
        elif 'INT16' == Type:
            Format = '<h'
        elif 'UINT32' == Type:
            Format = '<I'
        elif 'INT32' == Type:
            Format = '<i'
        elif 'UINT64' == Type:
            Format = '<Q'
        elif 'INT64' == Type:
            Format = '<q'
        elif 'FLOAT32' == Type:
            Format = '<f'
        elif 'FLOAT64' == Type:
            Format = '<d'
        elif 'REAL' == Type:
            Format = '<d'
        else:
            logging.error('ERROR26 Undefined Type!')

        return Format

    def defineFieldPositions(self):
        Blocks = self.info_Field.keys()

        for Block in Blocks:
            Fields = self.info_Field[Block].keys()
            index_Field = -1
            Field_previous = ''
            for Field in Fields:
                index_Field += 1
                if 0 == index_Field:
                    StartPosition = 0
                else:
                    StartPosition = self.positions_Field[Field_previous][1]
                Length = self.getFieldLength(Block, Field)
                self.positions_Field.update({Field: [StartPosition, StartPosition + Length]})
                Field_previous = Field

        return


    """
        followes a link chain and process the blocks on the links, step by step
    """
    def readLinks(self, nextentry, rblockid=None):
        result = []
        self.filepointer.seek(nextentry)
        self.actuallink = nextentry
        #reads the block_header and links from the actual block
        block_header, links = self.processHeaderAndLinkSections()
        while (nextentry != 0) and block_header is not None:
                #if the block is valid then process it
                blockvalue, nextentry = self.processBlocks(block_header, links, rblockid)
                del links
                del block_header
                result.append(blockvalue)
                if nextentry != 0:
                    #seek to the next block
                    self.filepointer.seek(nextentry)
                    self.actuallink = nextentry
                    block_header, links = self.processHeaderAndLinkSections()
                else:
                    block_header = None
        return result

    """
        Identificates a block's type and calls the proper block proessing method
        block_header, links is from processHeaderAndLinkSections()
        rblockid is a requested block type or type list
    """
    def processBlocks(self, block_header, links, rblockid=None):
        result = None
        nextentry = 0
        if block_header["id"][0:2] == "##":
            if rblockid is not None:
                #check for the requested block type
                if isinstance(rblockid, list):
                    ok = block_header["id"] in rblockid
                else:
                    ok = block_header["id"] == rblockid
            else:
                ok = True
            if ok:
                #if ok, lets find the proper block type
                if block_header["id"] == "##MD":
                    result, nextentry = self.blockMetaData(block_header, links)
                elif block_header["id"] == "##TX":
                    result, nextentry = self.blockTeXt(block_header, links)
                elif block_header["id"] == "##FH":
                    result, nextentry = self.blockFileHistory(block_header, links)
                elif block_header["id"] == "##CH":
                    result, nextentry = self.blockChannelHierarchy(block_header, links)
                elif block_header["id"] == "##AT":
                    result, nextentry = self.blockATtachment(block_header, links)
                elif block_header["id"] == "##EV":
                    result, nextentry = self.blockEVent(block_header, links)
                elif block_header["id"] == "##DG":
                    result, nextentry = self.blockDataGroup(block_header, links)
                elif block_header["id"] == "##CN":
                    result, nextentry = self.blockChaNnel(block_header, links)
                elif block_header["id"] == "##CG":
                    result, nextentry = self.blockChannelGroup(block_header, links)
                elif block_header["id"] == "##SI":
                    result, nextentry = self.blockSourceInfo(block_header, links)
                elif block_header["id"] == "##DT":
                    result, nextentry = self.blockDaTa(block_header, links)
                elif block_header["id"] == "##DV":
                    result, nextentry = self.blockDataValues(block_header, links)
                elif block_header["id"] == "##DL":
                    result, nextentry = self.blockDataList(block_header, links)
                elif block_header["id"] == "##LD":
                    result, nextentry = self.blockListData(block_header, links)
                elif block_header["id"] == "##HL":
                    result, nextentry = self.blockHeaderList(block_header, links)
            else:
                logging.error("ERROR5 Block id (%s) mismatch at %s, id is %s", str(rblockid), hex(self.actuallink), block_header["id"])
        else:
            logging.error("ERROR4 Invalid block at %s, id is %s", hex(self.actuallink), block_header["id"])
        return result, nextentry

    """
        Every mf4 file starts with some identification and helper fields
    """
    def blockID(self):
        # ASAMMDF documentation "4.5.1 Block Structure of IDBLOCK"
        self.filepointer.seek(0)
        data_bytes = self.filepointer.read(self.length_IDBLOCK)
        result = None
        id_file = self.readField(data_bytes, 'id_file')
        if id_file == 'MDF':
            result = {
                'id_file': id_file,
                'id_vers': self.readField(data_bytes, 'id_vers'),
                'id_prog': self.readField(data_bytes, 'id_prog'),
                'id_ver': self.readField(data_bytes, 'id_ver'),
                'id_unfin_flags': self.readField(data_bytes, 'id_unfin_flags'),
                'id_custom_unfin_flags': self.readField(data_bytes, 'id_custom_unfin_flags')
            }
        else:
            logging.error('ERROR-mf4reader-01 The file is not MDF!')

        if 400 > result['id_ver']:
            logging.error('ERROR-mf4reader-02 The MDF version is below 4.00!')

        return result

    """
        Reads the Header block from the mf4 file. This is the second block. It contains info about the measurement
        Called by the open()
    """
    def blockHeaDer(self):
        # ASAMMDF documentation "4.6 The Header Block HDBLOCK"
        link_HDBLOCK = self.length_IDBLOCK
        block_header, block_links_raw = self.readHeaderAndLinksSection(link_HDBLOCK)
        block_data = self.readDataSection(link_HDBLOCK, block_header)
        block_links = self.processLinks(block_header, block_links_raw, 0)

        result = {
            'header': block_header,
            'links': block_links,
            'data': block_data
        }

        return result

    """
        The MDBLOCK contains information encoded as XML string. For example, this can be
        comments for the measured data file, file history information or the identification of a channel.
        See documentation.
    """
    def blockMetaData(self, block_header, links):
        result = None
        data = self.content(block_header)
        if data is not None:
            #XML
            #TODO: parse xml
            result = data.decode("UTF8").strip("\0")
        return result, 0

    """
        The TXBLOCK is very similar to the MDBLOCK but only contains a plain string encoded in
        UTF-8. The text length results from the block size.
        See documentation.
    """
    def blockTeXt(self, block_header, links):
        # ASAMMDF documentation "4.8 The Text Block TXBLOCK"
        result = None
        data = self.content(block_header)
        if data is not None:
            result = data.decode('UTF8').strip('\0')
        return result, 0

    """
        The data section of the DTBLOCK contains a sequence of records. It contains records of all
        channel groups assigned to its parent DGBLOCK.
        See documentation.

        It is not readed by the initial file mapping. You can reach the data using the datagroup.getData()
    """
    def blockDaTa(self, block_header, links):
        # ASAMMDF documentation "4.21 The Data Block DTBLOCK"
        result = None
        data = self.content(block_header)
        if data is not None:
            result = {
                'header': block_header,
                'data':data
            }
        return result, 0

    """
        The data section of the DVBLOCK contains a sequence of records without their invalidation
        information. The size of these records is described by cg_data_bytes. It contains records of
        the (single) channel group assigned to its parent DGBLOCK. Unsorted writing is not allowed
        with this block. It is also not allowed to include a record id.
        See documentation.

        It is not readed by the initial file mapping. You can reach the data using the datagroup.getData()
    """
    def blockDataValues(self, block_header, links):
        result = None
        data = self.content(block_header)
        if data is not None:
            result = {
                    "datatype":"##DV",
                    "data":data
                    }
        return result, 0

    """
        The data section of the DIBLOCK contains a sequence of record invalidation information
        without the signal data. The size of these records is described by cg_inval_bytes. It contains
        records of the (single) channel group assigned to its parent DGBLOCK. Unsorted writing is
        not allowed with this block.
        See documentation.

        It is not readed by the initial file mapping. You can reach the data using the datagroup.getData()
    """
    def blockDataList(self, block_header, links):
        result = None
        links_processed = self.processLinks(block_header, links, 0)
        if 0 != links_processed['dl_dl_next']:
            logging.warning('The dl_dl_next != 0 is not handled yet!')

        dl_data = links_processed['dl_data']
        data_bytes = bytes()
        for dl_data_current in dl_data:
            self.setFilePosition(dl_data_current)
            header, links = self.processHeaderAndLinkSections()
            block_id = header['id']
            if '##DT' == block_id:
                result_blockDaTa, null = self.blockDaTa(header, links)
                data_bytes_current = result_blockDaTa['data']
                data_bytes += data_bytes_current
            else:
                logging.error('The ' + block_id + ' type block in DataList is not handled yet!')

        #data_bytes = self.content(block_header)
        #if data_bytes is not None:
        #    #TODO:link
        #    #datalist = self.blockDataList()
        result = {
            'datatype':block_header['id'],
            'data':data_bytes
        }

        return result, 0

    """
        The DLBLOCK references a list of data blocks (DTBLOCK) or a list of signal data blocks
        (SDBLOCK) or a list of reduction data blocks (RDBLOCK). This list of blocks is equivalent to
        using a single (signal/reduction) data block and can be used to avoid a huge data block by
        splitting it into smaller parts.
        See documentation.

        It is not readed by the initial file mapping. You can reach the data using the datagroup.getData()
    """
    def blockListData(self, block_header, links):
        result = None
        data = self.content(block_header)
        if data is not None:
            #TODO:link
            #datalist = self.blockDataList()
            result = {
                    "datatype":"##LD",
                    "data":""
                    }
        return result, 0

    """
        The HLBLOCK represents the "header" of a list of data blocks, i.e. the start of a linked list of
        DLBLOCKs. It contains information about the DLBLOCKs and the contained data blocks.
        See documentation.

        It is not readed by the initial file mapping. You can reach the data using the datagroup.getData()
    """
    def blockHeaderList(self, block_header, links):
        result = None
        data = self.content(block_header)
        if data is not None:
            #TODO:link
            #datalist = self.blockDataList()
            result = {
                    "datatype":"##HL",
                    "data":""
                    }
        return result, 0

    """
        The FHBLOCK describes who/which tool generated or changed the MDF file. Each
        FHBLOCK contains a change log entry for the MDF file.
        See documentation.
    """
    def blockFileHistory(self, block_header, links):
        # ASAMMDF documentation "4.9 The File History Block FHBLOCK"
        result = None
        nextentry = 0
        data = self.content(block_header)
        if data is not None:
            result = {
                'fh_time_ns':self.readField(data, 'fh_time_ns'),
                'fh_tz_offset_min':self.readField(data, 'fh_tz_offset_min'),
                'fh_dst_offset_min':self.readField(data, 'fh_dst_offset_min'),
                'fh_time_flags':self.readField(data, 'fh_time_flags')
            }

            links_processed = self.processLinks(block_header, links, 0)
            nextentry = links_processed['fh_fh_next']
        return result, nextentry

    """
        The ATBLOCK specifies attached data, either by referencing an external file or by
        embedding the data in the MDF file. Embedding the data, the data stream optionally can be
        compressed using the Deflate zip algorithm (see [12] and [24]).
        The ATBLOCK also holds information about the content of the attached data using a MIME
        content-type string (e.g. "video/avi" or "application/text", see [17] and [18]). For rules how to
        define custom MIME types please use the rules defined in ASAM ODS Version 5.1.1,
        chapter 8. In MDF, only the long form of the ASAM ODS MIME type string may be used, e.g.
        "application/x-asam.aolog".
        See documentation.
    """
    def blockATtachment(self, block_header, links):
        result = None
        nextentry = 0
        data = self.content(block_header, 36)
        if data is not None:
            result = {
                    "filename":self.getText(links[1]),
                    "mime":self.getText(links[2]),
                    "comment":self.getComment(links[3]),
                    #0.bit embedded data, 1.bit zip compressed embedded data, 2.bit valid checksum
                    "at_flags":DataType.readIntLittle(data, [0, DataType.INT16]),
                    "at_creator_index":DataType.readIntLittle(data, [2, DataType.INT16]),
                    #reserved 4 bytes
                    "at_md5_checksum":DataType.readChecksum(data, [4, DataType.INT8*16]),
                    "at_original_size":DataType.readIntLittle(data, [20, DataType.INT64]),
                    "at_embedded_size":DataType.readIntLittle(data, [28, DataType.INT64])
                    }
            nextentry = links[0]
        return result, nextentry

    """
        The EVBLOCK serves to describe an event. Each EVBLOCK stores a synchronization value
        to specify when the event occurred. Usually this will be a time stamp, but the event also can
        be synchronized with some other master channel type or the record index of a channel group
        (see 4.4.6 Synchronization Domains). Generally, an event defines a "point" in time or some
        other synchronization domain. For some event types, two "points" can be used to define a
        "range".
        See documentation.
    """
    def blockEVent(self, block_header, links):
        result = None
        nextentry = 0
        data = self.content(block_header)
        if data is not None:
            nextentry = links[0]
        return result, nextentry

    """
        The CHBLOCKs describe a logical ordering of the channels in a tree-like structure. This only
        serves to structure the channels and is totally independent to the data group and channel
        group structuring. A channel even may not be referenced at all or more than one time.
        See documentation.
    """
    def blockChannelHierarchy(self, FilePosition):
        StartingFilePosition = self.getFilePosition()
        self.setFilePosition(FilePosition)

        block_header, block_links = self.processHeaderAndLinkSections()

        result = None
        link_next = 0
        data_bytes = self.content(block_header)
        if data_bytes is not None:
            data_processed = {
                'ch_element_count': self.readField(data_bytes, 'ch_element_count'),
                'ch_type': self.readField(data_bytes, 'ch_type')
            }

            links_processed = self.processLinks(block_header, block_links, data_processed)
            link_next = links_processed['link_next']
            link_child_first = links_processed['link_child_first']
            name = self.getText(links_processed['ch_tx_name'])
            comment = self.getComment(links_processed['ch_md_comment'])

            if 0 != link_child_first:
                ch_childs = self.processChannelHierarchy(link_child_first)

            result = {
                'header': block_header,
                'links_processed': links_processed,
                'data_processed': data_processed,
                'name': name,
                'comment': comment,
                'ch_childs': ch_childs
            }

        self.setFilePosition(StartingFilePosition)

        return result, link_next

    """
        The DGBLOCK gathers information and links related to its data block. Thus the branch in the
        tree of MDF blocks that is opened by the DGBLOCK contains all information necessary to
        understand and decode the data block referenced by the DGBLOCK.
        The DGBLOCK can contain several channel groups. In this case the data group (and thus
        the MDF file) is "unsorted". If there is only one channel group in the DGBLOCK, the data
        group is "sorted"
        See documentation.
    """
    def blockDataGroup(self, block_header, links):
        # ASAMMDF documentation "4.13 The Data Group Block DGBLOCK"
        result = None
        nextentry = 0
        data_bytes = self.content(block_header)
        data_processed = None
        if data_bytes is not None:
            data_processed = {
                'dg_rec_id_size':self.readField(data_bytes, 'dg_rec_id_size')
            }

            links_processed = self.processLinks(block_header, links, data_processed)
            if 0 != links_processed['dg_cg_first']:
                self.links_Channelgroups.append(links_processed['dg_cg_first'])

            result = DataGroup(self,
                               data_processed['dg_rec_id_size'],
                               self.readLinks(links_processed['dg_cg_first'], '##CG'),
                               links_processed['dg_data'],
                               block_header,
                               links_processed,
                               data_processed,
                               self.getComment(links_processed['dg_md_comment'])
                               )

            nextentry = links_processed['dg_dg_next']
            if 0 != nextentry:
                self.links_Datagroups.append(nextentry)


        return result, nextentry

    """
        The CGBLOCK contains a collection of channels which are stored in one record, i.e. which
        have equal sampling.
        See documentation.
    """
    def blockChannelGroup(self, block_header, links):
        # ASAMMDF documentation "4.14 The Channel Group Block CGBLOCK"
        result = None
        nextentry = 0
        data_bytes = self.content(block_header)
        data_processed = None
        if data_bytes is not None:
            data_processed = {
                'cg_record_id':self.readField(data_bytes, 'cg_record_id'),
                'cg_cycle_count':self.readField(data_bytes, 'cg_cycle_count'),
                'cg_flags':self.readField(data_bytes, 'cg_flags'),
                'cg_path_separator':self.readField(data_bytes, 'cg_path_separator'),
                'cg_data_bytes':self.readField(data_bytes, 'cg_data_bytes'),
                'cg_inval_bytes':self.readField(data_bytes, 'cg_inval_bytes')
            }

            links_processed = self.processLinks(block_header, links, data_processed)

            channels = self.readLinks(links_processed['cg_cn_first'], '##CN')
            channel_0_Name = channels[0]['name']
            if 't' != channel_0_Name.lower():
                logging.warning('WARNING8 The 1st channel in Channels is not the Time!')

            result = {
                'channels': channels
            }
            #TODO: other links
            nextentry = links_processed['cg_cg_next']
            if 0 != nextentry:
                self.links_Channelgroups.append(nextentry)

        result_additional = {
            'header': block_header,
            'links_processed': links_processed,
            'data_processed': data_processed
        }
        result.update(result_additional)

        return result, nextentry

    """
        The CNBLOCK describes a channel, i.e. it contains information about the recorded signal
        and how its signal values are stored in the MDF file.
        See documentation.

        channel = signal
    """
    def blockChaNnel(self, block_header, links):
        # ASAMMDF documentation "4.16 The Channel Block CNBLOCK"
        result = None
        nextentry = 0
        data_bytes = self.content(block_header)
        if data_bytes is not None:
            data_processed = {
                'cn_type': self.readField(data_bytes, 'cn_type'),
                'cn_sync_type': self.readField(data_bytes, 'cn_sync_type'),
                'cn_data_type': self.readField(data_bytes, 'cn_data_type'),
                'cn_bit_offset': self.readField(data_bytes, 'cn_bit_offset'),
                'cn_byte_offset': self.readField(data_bytes, 'cn_byte_offset'),
                'cn_bit_count': self.readField(data_bytes, 'cn_bit_count'),
                'cn_flags': self.readField(data_bytes, 'cn_flags'),
                'cn_inval_bit_pos': self.readField(data_bytes, 'cn_inval_bit_pos'),
                'cn_precision': self.readField(data_bytes, 'cn_precision'),
                'cn_attachment_count': self.readField(data_bytes, 'cn_attachment_count'),
                'cn_val_range_min': self.readField(data_bytes, 'cn_val_range_min'),
                'cn_val_range_max': self.readField(data_bytes, 'cn_val_range_max'),
                'cn_limit_min': self.readField(data_bytes, 'cn_limit_min'),
                'cn_limit_max': self.readField(data_bytes, 'cn_limit_max'),
                'cn_limit_ext_min': self.readField(data_bytes, 'cn_limit_ext_min'),
                'cn_limit_ext_max': self.readField(data_bytes, 'cn_limit_ext_max')
            }

            links_processed = self.processLinks(block_header, links, data_processed)
            channel_name = self.getText(links_processed['cn_tx_name'])
            channel_comment = self.getComment(links_processed['cn_md_comment'])

            nextentry = links_processed['cn_cn_next']
            if 0 != nextentry:
                self.links_Channels.append(nextentry)

            if 0 != links_processed['cn_composition']:
                filepos_current = self.filepointer.tell()
                link_cn_composition = links_processed['cn_composition']
                self.filepointer.seek(link_cn_composition)
                block_header_cn_composition, links_cn_composition = self.processHeaderAndLinkSections()
                if '##CN' == block_header_cn_composition['id']:
                    result, null = self.blockChaNnel(block_header_cn_composition, links_cn_composition)
                    data_processed = result['data_processed']
                    links_processed = result['links_processed']
                    channel_name = self.getText(links_processed['cn_tx_name'])
                else:
                    logging.warning('WARNING7 The CABLOCK type cn_composition is not properly handled yet!')
                self.filepointer.seek(filepos_current)

            if None == result:
                result = {
                    'name':channel_name,
                    'source':self.getSource(links_processed['cn_si_source'])
                }
                #TODO: attachments, x axis links

                if 0 != links_processed['cn_cc_conversion']:
                    filepos_current = self.filepointer.tell()
                    link_cn_cc_conversion = links_processed['cn_cc_conversion']
                    self.filepointer.seek(link_cn_cc_conversion)
                    block_header_cn_cc_conversion, links_cn_cc_conversion = self.processHeaderAndLinkSections()
                    result_cn_cc_conversion, null = self.blockChannelConversion(block_header_cn_cc_conversion, links_cn_cc_conversion)
                    self.filepointer.seek(filepos_current)
                    result.update({'result_cn_cc_conversion': result_cn_cc_conversion})

                if 0 != links_processed['cn_si_source']:
                    filepos_current = self.filepointer.tell()
                    link_cn_si_source = links_processed['cn_si_source']
                    self.filepointer.seek(link_cn_si_source)
                    block_header_cn_si_source, links_cn_si_source = self.processHeaderAndLinkSections()
                    result_cn_si_source, null = self.blockSourceInfo(block_header_cn_si_source, links_cn_si_source)
                    self.filepointer.seek(filepos_current)
                    result.update({'result_cn_si_source': result_cn_si_source})


                result_additional = {
                    'header': block_header,
                    'links_processed': links_processed,
                    'data_processed': data_processed
                }
                result.update(result_additional)

        return result, nextentry

    """
        The CCBLOCK serves to specify a conversion formula to convert the raw values to
        physical values with a physical unit. The result of a conversion always is either
        a floating-point value (REAL) or a character string (UTF-8).
        See documentation.
    """
    def blockChannelConversion(self, block_header, links):
        # ASAMMDF documentation "4.17 The Channel Conversion Block CCBLOCK"
        result = None
        nextentry = 0
        data_bytes = self.content(block_header)
        if data_bytes is not None:
            data_processed = {
                'cc_type':self.readField(data_bytes, 'cc_type'),
                'cc_precision':self.readField(data_bytes, 'cc_precision'),
                'cc_flags':self.readField(data_bytes, 'cc_flags'),
                'cc_ref_count':self.readField(data_bytes, 'cc_ref_count'),
                'cc_val_count':self.readField(data_bytes, 'cc_val_count'),
                'cc_phy_range_min':self.readField(data_bytes, 'cc_phy_range_min'),
                'cc_phy_range_max':self.readField(data_bytes, 'cc_phy_range_max')
            }

            links_processed = self.processLinks(block_header, links, data_processed)

            cc_ref_count = data_processed['cc_ref_count']
            if 0 < cc_ref_count:
                cc_ref = links_processed['cc_ref']
                IsCCHeaderProcessed = False
                for cc_ref_current in cc_ref:
                    if 0 != cc_ref_current:
                        fileposition_current = self.getFilePosition()
                        self.setFilePosition(cc_ref_current)
                        cc_ref_header, cc_ref_links = self.processHeaderAndLinkSections()
                        if '##CC' == cc_ref_header['id']:
                            if IsCCHeaderProcessed:
                                logging.warning('WARNING10 Multiple CCBLOCKs were referenced in a CCBLOCK!')
                            IsCCHeaderProcessed = True
                            cc_ref_result, null = self.blockChannelConversion(cc_ref_header, cc_ref_links)
                            result = cc_ref_result
                        elif '##TX' == cc_ref_header['id']:
                            cc_ref_text, null = self.blockTeXt(cc_ref_header, cc_ref_links)
                            x = 1
                        self.setFilePosition(fileposition_current)

            cc_val_count = data_processed['cc_val_count']
            if 0 < cc_val_count:
                cc_val = []
                PositionInfo = self.positions_Field['cc_val']
                FieldLength = PositionInfo[1] - PositionInfo[0]
                for index_cc_val in range(0,cc_val_count):
                    StartPosition = PositionInfo[0] + index_cc_val * FieldLength
                    data_bytes_segment = data_bytes[StartPosition : (StartPosition + FieldLength)]
                    cc_val_current = self.readREAL(data_bytes_segment)
                    cc_val.append(cc_val_current)
                data_processed.update({'cc_val': cc_val})

            if None == result:
                result = {
                    'links_processed': links_processed,
                    'data_processed': data_processed
                }

        return result, nextentry

    """
        The SIBLOCK describes the source of an acquisition mode or of a signal. The source
        information is also used to ensure a unique identification of a channel
        See documentation.
    """
    def blockSourceInfo(self, block_header, links):
        # ASAMMDF documentation "4.15 The Source Information Block SIBLOCK"
        result = None
        nextentry = 0
        data_bytes = self.content(block_header)
        if data_bytes is not None:
            links_processed = self.processLinks(block_header, links, 0)

            data_processed = {
                'si_type': self.readField(data_bytes, 'si_type'),
                'si_bus_type': self.readField(data_bytes, 'si_bus_type'),
                'si_flags': self.readField(data_bytes, 'si_flags')
            }

            result = {
                'si_tx_name':self.getText(links_processed['si_tx_name']),
                'si_tx_path':self.getText(links_processed['si_tx_path']),
                'si_md_comment':self.getComment(links_processed['si_md_comment']),
                'si_type': data_processed['si_type'],
                'si_bus_type': data_processed['si_bus_type'],
                'si_flags': data_processed['si_flags']
            }

        return result, nextentry

    """
        reads the raw content of the block, if length is not defined, calculates the total datalength
    """
    def content(self, block_header, length=None):
        if length is None:
            length = block_header['length'] - (mf4reader.length_HeaderSection + block_header['link_count'] * DataType.LINK)
        if self.filesize >= self.filepointer.tell() + length:
            result = self.filepointer.read(length)
        else:
            logging.error('ERROR6 file overrun at %s', hex(self.filepointer.tell()))
            result = None
        return result

    """
        Processing the Header section and the links of the block
    """
    def processHeaderAndLinkSections(self):
        header = None
        links_raw = []
        # Checks can it read more mf4reader.length_HeaderSection bytes
        if self.filesize >= self.filepointer.tell() + mf4reader.length_HeaderSection:
            # Reads the actual block's header
            data = self.filepointer.read(mf4reader.length_HeaderSection)
            header = {
                "id":DataType.readStr(data, [0, 4*DataType.CHAR]), # Defines the type of block
                "length":DataType.readIntLittle(data, [8, 1*DataType.UINT64]), # Defines the full length of block in bytes
                "link_count":DataType.readIntLittle(data, [16, 1*DataType.UINT64]) # Defines the number of links_raw in the block
            }

            # Reads the associated links_raw of the block
            if header["link_count"] > 0:
                if self.filesize >= self.filepointer.tell() + header["link_count"] * DataType.LINK:
                    data = self.filepointer.read(header["link_count"] * DataType.LINK)
                    #a link is a byteposition of a block in the file
                    for i in range(0, header["link_count"]):
                        links_raw.append(DataType.readIntLittle(data, [i*DataType.LINK, DataType.LINK]))
                else:
                    logging.error("ERROR8 file overrun at %s", hex(self.filepointer.tell()))
        else:
            logging.error("ERROR7 file overrun at %s", hex(self.filepointer.tell()))
        return header, links_raw

    """
        Processing the Link section of the block
    """
    def processLinks(self, header, links_raw, data):
        # Processing/interpreting the links

        # Defining what the type of the Block
        if '##CC' == header['id']:
            # ASAMMDF documentation "4.17.1 Block Structure of CCBLOCK"
            cc_ref_count = data['cc_ref_count']
            cc_ref = []
            for index_cc_ref in range(0,cc_ref_count):
                cc_ref.append(links_raw[index_cc_ref + 4])

            links_processed = {
                'cc_tx_name': links_raw[0],
                'cc_md_unit': links_raw[1],
                'cc_md_comment': links_raw[2],
                'cc_cc_inverse': links_raw[3],
                'cc_ref': cc_ref
            }
            if (4 + cc_ref_count) < len(links_raw):
                logging.error('ERROR9 Something is not handled for CCBLOCK')

        elif '##CG' == header['id']:
            # ASAMMDF documentation "4.14.1 Block Structure of CGBLOCK"
            links_processed = {
                'cg_cg_next':links_raw[0],
                'cg_cn_first':links_raw[1],
                'cg_tx_acq_name':links_raw[2],
                'cg_si_acq_source':links_raw[3],
                'cg_sr_first':links_raw[4],
                'cg_md_comment':links_raw[5]
            }
            if 7 == len(links_raw):
                links_processed.update({'cg_cg_master': links_raw[6]})
                logging.warning('WARNING3 The handling of "remote master" is not yet implemented')

        elif '##CH' == header['id']:
            # ASAMMDF documentation "4.10.1 Block Structure of CHBLOCK"
            ch_element_count = data['ch_element_count']
            ch_element = []
            for index_ch_element in range(0,ch_element_count):
                link_DGBLOCK = links_raw[4 + (3*index_ch_element+1)]
                link_CGBLOCK = links_raw[4 + (3*index_ch_element+2)]
                link_CNBLOCK = links_raw[4 + (3*index_ch_element+3)]
                ch_element_current = {
                    'link_DGBLOCK': link_DGBLOCK,
                    'link_CGBLOCK': link_CGBLOCK,
                    'link_CNBLOCK': link_CNBLOCK
                }
                ch_element.append(ch_element_current)

            links_processed = {
                'link_next': links_raw[0],
                'link_child_first': links_raw[1],
                'ch_tx_name': links_raw[2],
                'ch_md_comment': links_raw[3],
                'ch_element': ch_element
            }


        elif '##CN' == header['id']:
            # ASAMMDF documentation "4.16.1 Block Structure of CNBLOCK"

            links_processed = {
                'cn_cn_next':links_raw[0],
                'cn_composition':links_raw[1],
                'cn_tx_name':links_raw[2],
                'cn_si_source':links_raw[3],
                'cn_cc_conversion':links_raw[4],
                'cn_data':links_raw[5],
                'cn_md_unit':links_raw[6],
                'cn_md_comment':links_raw[7]
            }

            if 0 != data:
                cn_attachment_count = data['cn_attachment_count']
                if (8 + cn_attachment_count) != len(links_raw):
                    logging.warning('WARNING4 The handling of "default X" flag is not yet implemented')
                if 0 < cn_attachment_count:
                    cn_at_reference = []
                    for i_cn_attachment_count in range(1,cn_attachment_count):
                        cn_at_reference_current = links_raw[7+i_cn_attachment_count]
                        cn_at_reference.append(cn_at_reference_current)

                    links_processed.update({'cn_at_reference': cn_at_reference})

        elif '##DG' == header['id']:
            # ASAMMDF documentation "4.13.1 Block Structure of DGBLOCK"
            links_processed = {
                'dg_dg_next':links_raw[0],
                'dg_cg_first':links_raw[1],
                'dg_data':links_raw[2],
                'dg_md_comment':links_raw[3]
            }

        elif '##DL' == header['id']:
            # ASAMMDF documentation "4.29.1 Block structure of DLBLOCK"
            dl_data = []
            link_count = header['link_count']
            for index_links in range(1,link_count):
                dl_data.append(links_raw[index_links])

            links_processed = {
                'dl_dl_next': links_raw[0],
                'dl_data': dl_data
            }


        elif '##FH' == header['id']:
            # ASAMMDF documentation "4.9.1 Block Structure of FHBLOCK"
            links_processed = {
                'fh_fh_next':links_raw[0],
                'fh_md_comment':links_raw[1]
            }

        elif '##HD' == header['id']:
            # ASAMMDF documentation "4.6.1 Block Structure of HDBLOCK"
            links_processed = {
                'hd_dg_first':links_raw[0],
                'hd_fh_first':links_raw[1],
                'hd_ch_first':links_raw[2],
                'hd_at_first':links_raw[3],
                'hd_ev_first':links_raw[4],
                'hd_md_comment':links_raw[5]
            }

        elif '##SI' == header['id']:
            # ASAMMDF documentation "4.15.1 Block Structure of SIBLOCK"
            links_processed = {
                'si_tx_name':links_raw[0],
                'si_tx_path':links_raw[1],
                'si_md_comment':links_raw[2]
            }

        elif '##TX' == header['id']:
            links_processed = {}

        else:
            logging.warning('WARNING09 The proper link processing of the ' + header['id'] + ' block is not yet implemented')
            links_processed = {}

        return links_processed

    def processChannelHierarchy(self, link_ch_first):
        ChannelHierarchy = []
        link_next = link_ch_first
        while 0 != link_next:
            ch_child, link_next = self.blockChannelHierarchy(link_next)
            ChannelHierarchy.append(ch_child)

        return ChannelHierarchy

    def mapDataGroups(self, link_dg_first):
        link_dg_next = link_dg_first
        #time_start = time.time()
        while 0 != link_dg_next:
            link_dg_next = self.readDataGroup(link_dg_next)
        #time_end = time.time()
        #time_elapsed_s = time_end - time_start
        #line_ElapsedTime = 'Elapsed time for Mapping: ' + str(time_elapsed_s) + ' s'
        #print(line_ElapsedTime)

        #lines_ChannelsList = []
        #lines_ChannelsList.append(line_ElapsedTime + '\n\n')
        #for Channel in self.ChannelsList:
        #    #if 1 < self.ChannelsList.count(Channel):
        #    #    logging.warning('WARNING19 This channel is several times in this list: ' + Channel)
        #    lines_ChannelsList.append(Channel + '\n')
        #file_ChannelsList = open((self.filepath[:-4] + '-ChannelsList_mf4reader.txt'), 'w')
        #file_ChannelsList.writelines(lines_ChannelsList)
        #file_ChannelsList.close()

        if TestMode:
            for Channel in self.ChannelsList:
                if 1 < self.ChannelsList.count(Channel):
                    print('This channel is in the measurement on that name multiple times: ' + Channel)

        return

    def readDataGroup(self, link_current):
        self.link_DGBLOCK_current = link_current

        self.filepointer.seek(link_current + self.length_HeaderSection)
        link_dg_next = DataType.readIntLittle(self.filepointer.read(DataType.LINK), [0, DataType.LINK])
        link_cg_next = DataType.readIntLittle(self.filepointer.read(DataType.LINK), [0, DataType.LINK])

        while 0 != link_cg_next:
            link_cg_next = self.readChannelGroup(link_cg_next)

        return link_dg_next

    def readChannelGroup(self, link_current):
        self.link_CGBLOCK_current = link_current

        self.filepointer.seek(link_current + self.length_HeaderSection)
        link_cg_next = DataType.readIntLittle(self.filepointer.read(DataType.LINK), [0, DataType.LINK])
        link_cn_next = DataType.readIntLittle(self.filepointer.read(DataType.LINK), [0, DataType.LINK])

        while 0 != link_cn_next:
            link_cn_next = self.readChaNnel(link_cn_next)

        return link_cg_next

    def readChaNnel(self, link_current):
        self.link_CNBLOCK_current = link_current

        self.filepointer.seek(link_current + self.length_HeaderSection)
        link_cn_next = DataType.readIntLittle(self.filepointer.read(DataType.LINK), [0, DataType.LINK])
        link_composition = DataType.readIntLittle(self.filepointer.read(DataType.LINK), [0, DataType.LINK])
        link_name = DataType.readIntLittle(self.filepointer.read(DataType.LINK), [0, DataType.LINK])

        IsCompositionContained = 0
        if 0 != link_composition:
            IsCompositionContained = 1

        Name = self.readTeXt(link_name)
        if     't' != Name.lower() \
           and 'CAN_DataFrame' not in Name \
           and 'CAN_ErrorFrame' not in Name \
           and 'CAN_RemoteFrame' not in Name:
            #if Name in self.ChannelsList:
            #    logging.warning('WARNING08 This channel is already in the list: ' + Name)
            self.ChannelsList.append(Name)
            ChannelInfo = [self.link_DGBLOCK_current, self.link_CGBLOCK_current, self.link_CNBLOCK_current, IsCompositionContained]
            self.ChannelsInfo.append(ChannelInfo)

        return link_cn_next

    def readComposition(self, link_current, Signal):
        # Reading and processing the composition(s)

        # Reading all the channels in this composition (channel)
        Channels = []
        ChannelsInfo_local = {}
        link_next = link_current
        while 0 != link_next:
            link_current = link_next
            header_block, links_block_raw = self.readHeaderAndLinksSection(link_current)
            links_block = self.processLinks(header_block, links_block_raw, 0)
            link_next = links_block['cn_cn_next']
            link_composition = links_block['cn_composition']
            link_name = links_block['cn_tx_name']
            name_channel = self.readTeXt(link_name)
            ChannelInfo = {
                'link_block': link_current,
                'link_composition': link_composition
            }
            Channels.append(name_channel)
            ChannelsInfo_local.update({name_channel: ChannelInfo})

        # Checking if the Signal is in the Channels
        if Signal in Channels:
            # The Channel of this Signal is in the current Channels list, NO further composition processing is needed

            # Checking if the Signal is several times in the Channels list
            if 1 < Channels.count(Signal):
                # The Signal is several times in the Channels list, therefore we can't be sure if the processed channel is the correct one
                logging.warning('WARNING20 This signal is several times in the composition: ' + Signal)

            # Getting the Channel info of the Signal's Channel and set the link to the Channel
            ChannelInfo = ChannelsInfo_local[Signal]
            link_channel = ChannelInfo['link_block']

        else:
            # The Channel of this Signal is NOT in the current Channels list, further composition processing is needed

            # Finding the Channel which contains the proper composition for the Signal
            IsParentChannelFound = False
            for Channel in Channels:
                if Channel in Signal:
                    # The Channel what contains the proper composition for the Signal is found
                    IsParentChannelFound = True
                    break

            # Checking if during the search the proper parent channel for the Signal is found
            if IsParentChannelFound:
                # The corrent parent channel for the Signal is found

                # Getting the channel info and processing the composition
                ChannelInfo = ChannelsInfo_local[Channel]
                link_composition = ChannelInfo['link_composition']
                header_composition = self.readHeaderSection(link_composition)

                # Checking if the composition is a Channel Block
                id_composition = header_composition['id']
                if '##CN' == id_composition:
                    # The composition is a Channel Block therefore process its
                    link_channel = self.readComposition(link_composition, Signal)

                else:
                    # The composition is not a Channel Block, what is not handled yet
                    logging.error('ERROR11 This type (' + id_composition + ') of composition is not handled yet! (' + Signal + ')')

            else:
                # No parent channel was found in the current channel
                logging.error('No proper parent Channel was found for that composed signal: ' + Signal)

        return link_channel

    def readREAL(self, data_bytes):
        Format = self.getUnpackFormat('REAL')
        Real = struct.unpack(Format, data_bytes)[0]
        return Real

    def readTeXt(self, link_text):
        self.filepointer.seek(link_text + self.length_ToLengthInHeader)
        Length_bytes = self.filepointer.read(self.length_UINT64)
        Length = struct.unpack('<Q', Length_bytes)[0]

        self.filepointer.seek(link_text + self.length_HeaderSection)
        Text_bytes = self.filepointer.read(self.length_BYTE * (Length - self.length_HeaderSection))
        Text = Text_bytes.decode('UTF8').strip('\0')

        return Text

    def readHeaderSection(self, link_block):
        self.filepointer.seek(link_block)
        bytes_header = self.filepointer.read(self.length_HeaderSection)
        header = {
            'id': self.readField(bytes_header, 'id'), # Defines the type of block
            'length': self.readField(bytes_header, 'length'), # Defines the full length of block in bytes
            'link_count': self.readField(bytes_header, 'link_count') # Defines the number of links_raw in the block
        }

        return header

    def readHeaderAndLinksSection(self, link_block):
        header = None
        links_raw = []

        self.filepointer.seek(link_block)
        # Checks can it read more mf4reader.length_HeaderSection bytes
        if self.filesize >= self.filepointer.tell() + self.length_HeaderSection:
            # Reads the actual block's header
            header = self.readHeaderSection(link_block)

            # Reads the associated links_raw of the block
            if header['link_count'] > 0:
                if self.filesize >= self.filepointer.tell() + header['link_count'] * self.length_LINK:
                    links_bytes = self.filepointer.read(header['link_count'] * self.length_LINK)
                    # A link is a byteposition of a block in the file
                    for i in range(0, header['link_count']):
                        link_bytes = links_bytes[i*self.length_LINK : (i*self.length_LINK + self.length_LINK)]
                        link_raw = self.readField(link_bytes, 'LINK')
                        links_raw.append(link_raw)
                else:
                    logging.error("ERROR8 file overrun at %s", hex(self.filepointer.tell()))
        else:
            logging.error("ERROR7 file overrun at %s", hex(self.filepointer.tell()))

        return header, links_raw

    def readDataBytes(self, link_data, length_data):
        self.filepointer.seek(link_data)

        if self.filesize >= self.filepointer.tell() + length_data:
            data_bytes = self.filepointer.read(length_data)
        else:
            logging.error('ERROR6 file overrun at %s', hex(self.filepointer.tell()))
            data_bytes = None

        return data_bytes

    def readField(self, data_bytes, Field):
        if    'id' == Field\
           or 'reserved' == Field\
           or 'length' == Field\
           or 'link_count' == Field:
            Block = 'Header'
        elif 'LINK' == Field:
            Block = 'Link'
        else:
            BlockId = Field[0:2]
            Block = BlockId.upper() + 'BLOCK'

        Type = self.info_Field[Block][Field][0]
        Count = self.info_Field[Block][Field][1]

        data_bytes_segment = data_bytes[self.positions_Field[Field][0] : self.positions_Field[Field][1]]
        if 'CHAR' == Type:
            data_Field = data_bytes_segment.decode("UTF8").strip()
        elif 'BYTE' == Type:
            logging.warning('WARNING22 Reading BYTE type of data field should not happen, something is wrong!')
            data_Field = 0
        elif 'LINK' == Type:
            logging.warning('WARNING22 Reading LINK type of data field should not happen, something is wrong!')
            data_Field = 0
        else:
            if 1 == Count:
                Format = self.getUnpackFormat(Type)
                data_Field = struct.unpack(Format, data_bytes_segment)[0]
            else:
                logging.error('ERROR28')

        return data_Field

    def readDataSection(self, link_block, header_block):
        id_block = header_block['id']

        link_count = header_block['link_count']
        length_HeaderAndLinkSection = self.length_HeaderSection + self.length_LINK*link_count
        length_block = header_block['length']

        link_data = link_block + length_HeaderAndLinkSection
        length_data = length_block - length_HeaderAndLinkSection
        data_bytes = self.readDataBytes(link_data, length_data)

        if '##HD' == id_block:
            # ASAMMDF documentation "4.6.1 Block Structure of HDBLOCK"
            data = {
                'hd_start_time_ns':self.readField(data_bytes, 'hd_start_time_ns'),
                'hd_start_time_ns':self.readField(data_bytes, 'hd_start_time_ns'),
                'hd_tz_offset_min':self.readField(data_bytes, 'hd_tz_offset_min'),
                'hd_dst_offset_min':self.readField(data_bytes, 'hd_dst_offset_min'),
                'hd_time_flags':self.readField(data_bytes, 'hd_time_flags'),
                'hd_time_class':self.readField(data_bytes, 'hd_time_class'),
                'hd_flags':self.readField(data_bytes, 'hd_flags'),
                'hd_start_angle_rad':self.readField(data_bytes, 'hd_start_angle_rad'),
                'hd_start_distance_m':self.readField(data_bytes, 'hd_start_distance_m')
            }

        elif '##DG' == id_block:
            data = {
                'dg_rec_id_size':self.readField(data_bytes, 'dg_rec_id_size')
            }

        elif '##CG' == id_block:
            data = {
                'cg_record_id':self.readField(data_bytes, 'cg_record_id'),
                'cg_cycle_count':self.readField(data_bytes, 'cg_cycle_count'),
                'cg_flags':self.readField(data_bytes, 'cg_flags'),
                'cg_path_separator':self.readField(data_bytes, 'cg_path_separator'),
                'cg_data_bytes':self.readField(data_bytes, 'cg_data_bytes'),
                'cg_inval_bytes':self.readField(data_bytes, 'cg_inval_bytes')
            }

        elif '##CN' == id_block:
            data = {
                'cn_type':self.readField(data_bytes, 'cn_type'),
                'cn_sync_type':self.readField(data_bytes, 'cn_sync_type'),
                'cn_data_type':self.readField(data_bytes, 'cn_data_type'),
                'cn_bit_offset':self.readField(data_bytes, 'cn_bit_offset'),
                'cn_byte_offset':self.readField(data_bytes, 'cn_byte_offset'),
                'cn_bit_count':self.readField(data_bytes, 'cn_bit_count'),
                'cn_flags':self.readField(data_bytes, 'cn_flags'),
                'cn_inval_bit_pos':self.readField(data_bytes, 'cn_inval_bit_pos'),
                'cn_precision':self.readField(data_bytes, 'cn_precision'),
                'cn_attachment_count':self.readField(data_bytes, 'cn_attachment_count'),
                'cn_val_range_min':self.readField(data_bytes, 'cn_val_range_min'),
                'cn_val_range_max':self.readField(data_bytes, 'cn_val_range_max'),
                'cn_limit_min':self.readField(data_bytes, 'cn_limit_min'),
                'cn_limit_max':self.readField(data_bytes, 'cn_limit_max'),
                'cn_limit_ext_min':self.readField(data_bytes, 'cn_limit_ext_min'),
                'cn_limit_ext_max':self.readField(data_bytes, 'cn_limit_ext_max')
            }

            cn_flags = data['cn_flags']
            if 0 != cn_flags:
                header_cn_block, links_cn_block_raw = self.readHeaderAndLinksSection(link_block)
                links_cn_block = self.processLinks(header_cn_block, links_cn_block_raw, data)
                channel_name = self.readTeXt(links_cn_block['cn_tx_name'])

                cn_flags_bin = bin(cn_flags)
                while 17 != len(cn_flags_bin):
                    cn_flags_bin = cn_flags_bin[:2] + '0' + cn_flags_bin[2:]
                if TestMode:
                    logging.info('CN Flags (' + cn_flags_bin + ') is set for ' + channel_name)

                if '1' == cn_flags_bin[-1]:
                    if TestMode:
                        logging.info('CN Invalid flag is set for ' + channel_name)

                if '1' == cn_flags_bin[-3]:
                    cn_precision_hex = hex(data['cn_precision'])
                    if TestMode:
                        logging.info('CN Precision is ' + cn_precision_hex + ' for ' + channel_name)

                if '1' == cn_flags_bin[-4]:
                    if TestMode:
                        logging.info('CN Min is ' + str(data['cn_val_range_min']) + ', CN Max is ' + str(data['cn_val_range_max']) + ' for ' + channel_name)

        elif '##CC' == id_block:
            data = {
                'cc_type':self.readField(data_bytes, 'cc_type'),
                'cc_precision':self.readField(data_bytes, 'cc_precision'),
                'cc_flags':self.readField(data_bytes, 'cc_flags'),
                'cc_ref_count':self.readField(data_bytes, 'cc_ref_count'),
                'cc_val_count':self.readField(data_bytes, 'cc_val_count'),
                'cc_phy_range_min':self.readField(data_bytes, 'cc_phy_range_min'),
                'cc_phy_range_max':self.readField(data_bytes, 'cc_phy_range_max')
            }

            cc_val_count = data['cc_val_count']
            if 0 < cc_val_count:
                cc_val = []
                PositionInfo = self.positions_Field['cc_val']
                FieldLength = PositionInfo[1] - PositionInfo[0]
                for index_cc_val in range(0,cc_val_count):
                    StartPosition = PositionInfo[0] + index_cc_val * FieldLength
                    data_bytes_segment = data_bytes[StartPosition : (StartPosition + FieldLength)]
                    cc_val_current = self.readREAL(data_bytes_segment)
                    cc_val.append(cc_val_current)
                data.update({'cc_val': cc_val})

            cc_flags = data['cc_flags']
            if 0 != cc_flags:
                cc_flags_bin = bin(cc_flags)
                while 5 != len(cc_flags_bin):
                    cc_flags_bin = cc_flags_bin[:2] + '0' + cc_flags_bin[2:]
                if TestMode:
                    logging.info('CC Flags (' + cc_flags_bin + ') is set for ' + self.Channel_current)

                if '1' == cc_flags_bin[4]:
                    if TestMode:
                        logging.info('CC Precision is ' + str(data['cc_precision']) + ' for ' + self.Channel_current)

                if '1' == cc_flags_bin[3]:
                    if TestMode:
                        logging.info('CC Min is ' + str(data['cc_phy_range_min']) + ', CC Max is ' + str(data['cc_phy_range_max']) + ' for ' + self.Channel_current)


        elif '##TX' == id_block:
            tx_data = self.readTeXt(link_block)
            data = {
                'tx_data': tx_data
            }

        else:
            logging.warning('WARNING17 This data processing for this block type (' + id_block + ') is not handled yet!')
            data = {}

        return data

    def readBlock(self, link_block):
        # Reading the Block info

        # Reading the Header (and ID), Links and Data of the Block
        header_block, links_block_raw = self.readHeaderAndLinksSection(link_block)
        id_block = header_block['id']
        data_block = self.readDataSection(link_block, header_block)
        links_block = self.processLinks(header_block, links_block_raw, data_block)

        # Creating the returning info of the Block
        info_block = {
            'id_block': id_block,
            'header_block': header_block,
            'links_block': links_block,
            'data_block': data_block
        }

        # Doing some special processing for specific Block types
        if '##CN' == id_block:
            # ASAMMDF documentation "4.16 The Channel Block CNBLOCK"

            channel_name = self.readTeXt(links_block['cn_tx_name'])
            self.Channel_current = channel_name

            # Getting the Conversion value for Channel Block

            # Getting the Conversion link and processing it
            link_conversion = links_block['cn_cc_conversion']
            if 0 != link_conversion:
                # The Conversion link is not 0, Conversion is present

                # Processing the Conversion
                info_conversion = self.readBlock(link_conversion)

                # Checking if the returning Conversion value (offset and factor) is in proper format
                cc_val = info_conversion['data_block']['cc_val']
                if 2 == len(cc_val):
                    # The returning Conversion value is in proper format
                    conversion = {
                        'offset': cc_val[0],
                        'factor': cc_val[1]
                    }

                else:
                    # The returning Conversion value is NOT in proper format, setting it to default
                    conversion = {
                        'offset': 0,
                        'factor': 1
                    }

                    logging.warning('WARNING03 Something is wrong with the conversion of this channel, it was set to default: ' + channel_name)

            else:
                # The Conversion link is 0, Conversion is set to default
                conversion = {
                    'offset': 0,
                    'factor': 1
                }

            info_block.update({'conversion': conversion})

        elif '##CC' == id_block:
            # ASAMMDF documentation "4.17 The Channel Conversion Block CCBLOCK"
            # Processing the Conversion value of a Conversion Block

            # Getting the number of referenced conversions and the type of the Conversion
            cc_ref_count = data_block['cc_ref_count']
            cc_type = data_block['cc_type']

            # Checking the type of the Conversion
            if 0 == cc_type:
                # Type 0 = 1:1 conversion -> default conversion
                cc_val = [0, 1]
                data_block.update({'cc_val': cc_val})

            elif 1 == cc_type:
                # Type 1 = linear conversion

                # Checking if the Conversion value is in proper format
                cc_val = data_block['cc_val']
                if 2 != len(cc_val):
                    # The conversion value is not in proper format, something is wrong!
                    logging.warning('WARNING05 The CC value format for conversion type 1 is not proper!')

            elif 7 == cc_type:
                # Type 7 = value to text/scale conversion tabular look-up

                # Checking the number of references
                if 0 < cc_ref_count:
                    # Processing the Conversion references
                    links_cc_ref = links_block['cc_ref']
                    IsCCHeaderProcessed = False
                    for link_cc_ref in links_cc_ref:
                        if 0 != link_cc_ref:
                            # Processing the Conversion reference
                            info_cc_ref = self.readBlock(link_cc_ref)
                            id_cc_ref = info_cc_ref['id_block']
                            if '##CC' == id_cc_ref:
                                # The referenced link is CCBLOCK therefore its value should be used
                                if IsCCHeaderProcessed:
                                    logging.warning('WARNING01 Multiple CCBLOCKs were referenced in a Conversion!')
                                IsCCHeaderProcessed = True
                                info_block = info_cc_ref
                            elif '##TX' == id_cc_ref:
                                # The referenced link is a TXBLOCK what is not handled yet
                                pass # logging.warning('WARNING02 Enum texts are not yet handled in a CCBLOCK!')

                    if False == IsCCHeaderProcessed:
                        logging.warning('No value was referenced for the current conversion')

                else:
                    logging.warning('The number of references for Type 7 Conversion is 0 what should not happen!')
            else:
                logging.warning('WARNING04 This conversion type (' + str(cc_type) + ') is not handled yet!')

        return info_block

    """
        closes the mf4 file, must call if not used anymore
    """
    def close(self):
        self.filepointer.close()

    """
        *** These methods used by other methods to get a linked field's value
        *** BEGIN
    """

    #higher level block query
    def getBlock(self, link, rblockid):
        result = None
        if link != 0:
            self.filepointer.seek(link)
            block_header, links = self.processHeaderAndLinkSections()
            if block_header is not None:
                blockvalue, nextentry = self.processBlocks(block_header, links, rblockid)
                result = blockvalue
        return result

    #Reads a comment block (text or metadata blocks)
    def getComment(self, link):
        blockvalue = self.getBlock(link, ["##TX", "##MD"])
        if blockvalue is not None:
            result = blockvalue
        else:
            result = ""
        return result

    #Reads a text block
    def getText(self, link):
        blockvalue = self.getBlock(link, "##TX")
        if blockvalue is not None:
            result = blockvalue
        else:
            result = ""
        return result

    #Reads a metadata block
    def getMeta(self, link):
        blockvalue = self.getBlock(link, "##MD")
        if blockvalue is not None:
            result = blockvalue
        else:
            result = ""
        return result

    #Reads a sourceinfo block
    def getSource(self, link):
        return self.getBlock(link, "##SI")

    # Give back the StartPosition and the EndPosition
    def getPositions(self, Field):
        StartPosition = self.StartPosition[Field]
        EndPosition = StartPosition + mf4reader.length_Field[Field]
        return [StartPosition, EndPosition]

    # Give back the current file position
    def getFilePosition(self):
        return self.filepointer.tell()

    # Set the the current file position
    def setFilePosition(self, position):
        self.filepointer.seek(position)
        return

    def readDTBLOCK(self, link_data, info_rawData):
        # Reading a simple Data block
        # Getting the raw data for the Time and Signal channel

        # Defining the length of the Data section
        header_block = self.readHeaderSection(link_data)
        length_data = header_block['length'] - self.length_HeaderSection

        # Checking if the bytes length of the data segment is fitting in the Data section
        cg_data_bytes = info_rawData['cg_data_bytes']
        if 0 == length_data % cg_data_bytes:
            # The bytes length of the data segment is fitting in the Data section so process the data

            # Getting the Time and Signal info:
            #   - byte length of the Time (cn_byte_count_Time)
            #   - byte offset of the Value (cn_byte_offset_Value)
            #   - byte length of the Value (cn_byte_count_Value)
            cn_byte_count_Time = info_rawData['cn_byte_count_Time']
            cn_byte_offset_Value = info_rawData['cn_byte_offset_Value']
            cn_byte_count_Value = info_rawData['cn_byte_count_Value']

            # Processing the data and collecting the relevant data segments
            data_preprocessed = []
            # Defining the null point for the data pointer (where the data section starts in the file)
            pointer_data_null = link_data + self.length_HeaderSection
            # Defining the number of cycles in the Data section
            count_cycles = int(length_data / cg_data_bytes)
            for index in range(0,count_cycles):
                # Getting the raw data of Time
                self.filepointer.seek(pointer_data_null + (index * cg_data_bytes))
                data_bytes_Time = self.filepointer.read(cn_byte_count_Time)

                # Getting the raw data of Value
                self.filepointer.seek(pointer_data_null + ((index * cg_data_bytes) + cn_byte_offset_Value))
                data_bytes_Value = self.filepointer.read(cn_byte_count_Value)

                # Expanding the raw data array
                data_preprocessed.append([data_bytes_Time, data_bytes_Value])

        else:
            # The bytes length of the data segment is NOT fitting in the Data section, something is wrong!
            logging.error('ERROR23 The length of the DTBLOCK is not proper according to the data bytes length')

        return data_preprocessed

    def getRawData(self, link_data, info_rawData):
        # Getting the raw data with some preprocessing
        # The preprocessing is getting just those data segments what contain Time and Signal data

        # Getting the header info and links and checking whether it is simple Data block or Data List block
        header_block, links_block_raw = self.readHeaderAndLinksSection(link_data)
        block_id = header_block['id']
        if '##DT' == block_id:
            # It's a simple Data block so just get the raw data
            data_bytes = self.readDTBLOCK(link_data, info_rawData)

        elif '##DL' == block_id:
            # It's a Data List block so get the raw data from all referenced Data blocks and attach them together

            # Processing the links to interpret them
            links_processed = self.processLinks(header_block, links_block_raw, 0)

            # Checking if next Data List is defined, if yes, it's unknown for us yet
            if 0 != links_processed['dl_dl_next']:
                logging.warning('The dl_dl_next != 0 is not handled yet!')

            # Getting the Data blocks defined in this Data List and processing them
            links_dl_data = links_processed['dl_data']
            data_bytes = []
            for link_dl_data in links_dl_data:
                # Processing each reference Data block one-by-one

                # Getting the header info of the block and checking if it is simple Data block
                header_dl_block = self.readHeaderSection(link_dl_data)
                if '##DT' == header_dl_block['id']:
                    # It's a simple Data block, process it and expand the Data list
                    data_bytes_current = self.readDTBLOCK(link_dl_data, info_rawData)
                    data_bytes += data_bytes_current

                else:
                    # It's NOT a simple Data block, what is not handled yet!
                    logging.error('The ' + block_id + ' type block in DataList is not handled yet!')

        else:
            # It's either a simple Data block nor a Data List block and it's NOT handled yet
            logging.warning('The getting of raw data for '+ block_id + ' data block type is not yet implemented!')

        return data_bytes

    def getByteCount(self, bit_count):
        # Defining the byte length (byte_count) in the bit length (bit_count) fit in

        if 0 == bit_count % 8:
            byte_count = int(bit_count / 8)
        else:
            byte_count = int((bit_count + (8 - bit_count % 8)) / 8)

        return byte_count

    def getSignalValue(self, SignalName):
        # Checking if the Signal is in the processed channels list
        # If yes then the proper links for the data are directly available
        # If not then further processing is needed for getting the links
        if SignalName in self.ChannelsList:
            if 1 < self.ChannelsList.count(SignalName):
                index_Signal = -1
                for index in range(0, self.ChannelsList.count(SignalName)-1):
                    index_Signal = self.ChannelsList.index(SignalName, index_Signal+1)
                    ChannelInfo = self.ChannelsInfo[index_Signal]
                    link_DGBLOCK = ChannelInfo[0]
                    block_DGBLOCK = self.readBlock(link_DGBLOCK)
                    link_data = block_DGBLOCK['links_block']['dg_data']
                    if 0 != link_data:
                        break
            else:
                index_Signal = self.ChannelsList.index(SignalName)
                ChannelInfo = self.ChannelsInfo[index_Signal]
                link_DGBLOCK = ChannelInfo[0]
            link_CGBLOCK = ChannelInfo[1]
            link_CNBLOCK = ChannelInfo[2]

        else:
            # The Signal is not in the processed channels list therefore further processing is needed for getting the proper links
            # The Signal is probably in a channel's composition

            # Checking and counting the number of dots in the Signal. If dot is present it shows it's probably a composition
            count_dots = SignalName.count('.')
            if 0 != count_dots:
                # Dot is present in the Signal name therefore the Signal is probably a channel's composition

                # Finding the channel in the composition could be.
                # The name of the channel, where the Signal composition is, should be a part of the Signal's name,
                # the possible parts are separated with dots.
                # Checking all the possible parts (starting from the beginning and also from the ending, separated with dots)
                # what can be the channel in the compostion could be.
                IsMainTagInChannelsList = False

                # Checking all the possible parts from the beginning if it is in the channels list
                mainTag = ''
                sideTag = SignalName
                for index in range(0, sideTag.count('.')):
                    position_firstDot = sideTag.find('.')
                    mainTag = SignalName[:position_firstDot]
                    sideTag = SignalName[position_firstDot+1:]
                    if mainTag in self.ChannelsList:
                        # The part is found in the channels list, the Signal is a composition of the channel on that name
                        IsMainTagInChannelsList = True
                        break

                # Checking if the channel is found in the proper composition (Signal) is
                if False == IsMainTagInChannelsList:
                    # Checking all the possible parts from the ending if it is in the channels list
                    mainTag = ''
                    sideTag = SignalName
                    for index in range(0, sideTag.count('.')):
                        position_lastDot = sideTag.rfind('.')
                        mainTag = SignalName[position_lastDot+1:]
                        sideTag = SignalName[:position_lastDot]
                        if mainTag in self.ChannelsList:
                            # The part is found in the channels list, the Signal is a composition of the channel on that name
                            IsMainTagInChannelsList = True
                            break

                if False == IsMainTagInChannelsList:
                    # None of the possible parts is in the channels list, so the composition hierarchy is not proper
                    # or the Signal is not in the measurement
                    logging.error('ERROR20 This signal is NOT a composition: ' + SignalName)

                else:
                    # The channel where the Signal's composition is found

                    # Getting the channel info of the composition's parent
                    index_mainTag = self.ChannelsList.index(mainTag)
                    ChannelInfo = self.ChannelsInfo[index_mainTag]
                    link_DGBLOCK = ChannelInfo[0]
                    link_CGBLOCK = ChannelInfo[1]
                    link_CNBLOCK = ChannelInfo[2]

                    # Checking if this channel is containing composition
                    IsCompositionContained = ChannelInfo[3]
                    if 0 == IsCompositionContained:
                        # This channel is not containing composition, so the composition hierarchy is wrong
                        # or the Signal is not in the measurement on that name
                        logging.error('ERROR19 This signal is NOT a composition: ' + SignalName)

                    # Getting the link to the composition and processing it
                    header_CNBLOCK, links_CNBLOCK_raw = self.readHeaderAndLinksSection(link_CNBLOCK)
                    links_CNBLOCK = self.processLinks(header_CNBLOCK, links_CNBLOCK_raw, 0)
                    link_CNBLOCK_composition = links_CNBLOCK['cn_composition']

                    # Checking if the composition link is proper
                    if 0 == link_CNBLOCK_composition:
                        # The composition link is not proper, there was something error during the mapping or during the composition processing
                        logging.error('ERROR27 This signal should have composition but it is not!: ' + SignalName)

                    # Getting the link of the channel containing the Signal
                    link_CNBLOCK = self.readComposition(link_CNBLOCK_composition, SignalName)

            else:
                # There is no dot in the Signal name therefore the Signal is not a channel's composition
                # therefore the Signal is not in the measurement on this name
                logging.error('ERROR17 This signal is NOT a composition: ' + SignalName)

        # Reading the Datagroup, Channelgroup and Channel block of the Signal and the proper Channel block of the Time
        info_DGBLOCK = self.readBlock(link_DGBLOCK)
        info_CGBLOCK = self.readBlock(link_CGBLOCK)
        info_CNBLOCK = self.readBlock(link_CNBLOCK)
        info_CNBLOCK_time = self.readBlock(info_CGBLOCK['links_block']['cg_cn_first'])

        # Getting the Channelgroup info:
        #   - byte length of the data segment (cg_data_bytes)
        #   - the number of cycles (cg_cycle_count)
        cg_data_bytes = info_CGBLOCK['data_block']['cg_data_bytes']
        cg_cycle_count = info_CGBLOCK['data_block']['cg_cycle_count']

        # Getting the Time Channel info:
        #   - conversion offset (cc_val_offset_Time)
        #   - conversion factor (cc_val_factor_Time)
        #   - bit length of the signal (cn_bit_count_Time)
        #   - byte length of the signal, in the bit length is fit (cn_byte_count_Time)
        cc_val_offset_Time = info_CNBLOCK_time['conversion']['offset']
        cc_val_factor_Time = info_CNBLOCK_time['conversion']['factor']
        cn_bit_count_Time = info_CNBLOCK_time['data_block']['cn_bit_count']
        cn_byte_count_Time = self.getByteCount(cn_bit_count_Time)

        # Getting the Value Channel info:
        #   - conversion offset (cc_val_offset_Value)
        #   - conversion factor (cc_val_factor_Value)
        #   - data type (cn_data_type_Value)
        #   - byte offset in the data segment (cn_byte_offset_Value)
        #   - bit length of the signal (cn_bit_count_Time)
        #   - byte length of the signal, in the bit length is fit (cn_byte_count_Time)
        #   - bit offset in the data byte segment, if the bit length not fill a byte length (cn_bit_offset_Value)
        cc_val_offset_Value = info_CNBLOCK['conversion']['offset']
        cc_val_factor_Value = info_CNBLOCK['conversion']['factor']
        cn_data_type_Value = info_CNBLOCK['data_block']['cn_data_type']
        cn_byte_offset_Value = info_CNBLOCK['data_block']['cn_byte_offset']
        cn_bit_count_Value = info_CNBLOCK['data_block']['cn_bit_count']
        cn_byte_count_Value = self.getByteCount(cn_bit_count_Value)
        cn_bit_offset_Value = info_CNBLOCK['data_block']['cn_bit_offset']

        # Getting the raw data from the Datagroup
        link_data = info_DGBLOCK['links_block']['dg_data']
        #time_start = time.time()
        info_rawData = {
            'cg_data_bytes': cg_data_bytes,
            'cn_byte_count_Time': cn_byte_count_Time,
            'cn_byte_offset_Value': cn_byte_offset_Value,
            'cn_byte_count_Value': cn_byte_count_Value,
        }
        data_bytes = self.getRawData(link_data, info_rawData)
        #time_end = time.time()
        #time_elapsed_s = time_end - time_start
        #line_ElapsedTime = 'Elapsed time for Data read: ' + str(time_elapsed_s) + ' s'
        #print(line_ElapsedTime)

        values = []
        values_Time = []
        values_Value = []
        #lines_CSV = []
        #lines_CSV.append('T,'+ SignalName + '\n')
        #lines_CSV.append('s,\n')
        if cg_cycle_count == len(data_bytes):
            for i in range(0,cg_cycle_count):
                value_Time = DataType.readFloat64(data_bytes[i][0], [0, cn_byte_count_Time]) * cc_val_factor_Time + cc_val_offset_Time
                values_Time.append(value_Time)

                # ASAMMDF documentation "4.16.1 Block Structure of CNBLOCK", cn_data_type
                if 0 == cn_data_type_Value or 2 == cn_data_type_Value or 4 == cn_data_type_Value:
                    # LE Byte order is used
                    struct_unpack_format = '<'
                elif 1 == cn_data_type_Value or 3 == cn_data_type_Value or 5 == cn_data_type_Value:
                    # BE Byte order is used
                    struct_unpack_format = '>'
                else:
                    logging.error('ERROR11 The "cn_data_type = ' + cn_data_type_Value + '" signal data type reading is not handled yet! (Signal: ' + SignalName + ')')

                datasegment_Value = data_bytes[i][1]
                if 0 == cn_data_type_Value or 1 == cn_data_type_Value or 2 == cn_data_type_Value or 3 == cn_data_type_Value:
                    # Integer data type

                    if 1 == cn_byte_count_Value:
                        struct_unpack_format += 'B'
                    elif 2 == cn_byte_count_Value:
                        struct_unpack_format += 'H'
                    elif 4 == cn_byte_count_Value:
                        struct_unpack_format += 'L'
                    elif 8 == cn_byte_count_Value:
                        struct_unpack_format += 'Q'
                    else:
                        logging.error('ERROR12 The "cn_byte_count_Value = ' + str(cn_byte_count_Value) + '" signal length reading for integer is not handled yet! (Signal: ' + SignalName + ')')

                    if 2 == cn_data_type_Value or 3 == cn_data_type_Value:
                        # Signed data type
                        struct_unpack_format = struct_unpack_format.lower()

                    value_Value_raw = struct.unpack(struct_unpack_format, datasegment_Value)[0]
                    if 0 != cn_bit_count_Value % 8:
                        if 2 == cn_data_type_Value or 3 == cn_data_type_Value:
                            logging.warning('WARNING6 The value processing for signed integer non-byte-length signal is not tested yet! (Signal: ' + SignalName + ')')

                        value_Value_raw_bin = bin(value_Value_raw)[2:]
                        while (cn_byte_count_Value * 8) > len(value_Value_raw_bin):
                            value_Value_raw_bin = '0' + value_Value_raw_bin

                        if 0 == cn_bit_offset_Value:
                            value_Value_raw_bin_cut = value_Value_raw_bin[-(cn_bit_offset_Value+cn_bit_count_Value) : ]
                        else:
                            value_Value_raw_bin_cut = value_Value_raw_bin[-(cn_bit_offset_Value+cn_bit_count_Value) : -cn_bit_offset_Value]

                        value_Value_raw = int('0b' + value_Value_raw_bin_cut, 2)

                elif 4 == cn_data_type_Value or 5 == cn_data_type_Value:
                    # Floating-point data type
                    if (4 * 8) == cn_bit_count_Value:
                        struct_unpack_format += 'f'
                        cn_byte_count_Value = 4
                    elif (8 * 8) == cn_bit_count_Value:
                        struct_unpack_format += 'd'
                        cn_byte_count_Value = 8
                    else:
                        logging.error('ERROR13 The signal value processing "cn_bit_count = ' + cn_bit_count_Value + '" signal length for float is not handled yet! (Signal: ' + SignalName + ')')

                    value_Value_raw = struct.unpack(struct_unpack_format, datasegment_Value)[0]

                value_Value = value_Value_raw * cc_val_factor_Value + cc_val_offset_Value
                values_Value.append(value_Value)

                values.append([value_Time,value_Value])

                #value_Time_CSV = format(value_Time, '.6f')
                #value_Value_CSV = format(value_Value, '.6f')
                #line_CSV = str(value_Time_CSV) + ',' + str(value_Value_CSV) + '\n'
                #lines_CSV.append(line_CSV)

            #filepath = self.filepath
            #position_LastPer = filepath.rfind('/')
            #if -1 == position_LastPer:
            #    path_FileDirectory = ''
            #else:
            #    path_FileDirectory = filepath[:position_LastPer+1]
            #f_test_out = open(path_FileDirectory + 'out-' + SignalName + '.csv', 'w')
            #f_test_out.writelines(lines_CSV)
            #f_test_out.close()

        else:
            logging.error('ERROR24 The raw data reading and preprocessing was not proper! (' + SignalName + ')')

        Signal = {'Time': values_Time, 'Value': values_Value}
        return Signal

    """
        *** These methods used by other methods to get a linked field's value
        *** END
    """

if TestMode:
    MF4Test = mf4reader()
    MF4Test.open('libs/mf4extractor/test/PJHF_2020-12-25_14-44_40.MF4')
    MF4Test.getSignalValue('YawRateActMeas.ActManLat')
    print('MF4Test was executed!')
