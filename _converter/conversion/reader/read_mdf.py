"""
@author: jez2bp
"""

import re
import numpy
import logging
import asammdf
import _converter.messages.messages


class ReadMDF:
    """The ReadMDF class reads the MDF file and contructs an asamdf.MDF object
    which contains only the signals what was definied as plt_signals.

    -Init:
        - filepath: the absulute path of an MDF file
        - plt_signals: list of strings with the name of the needed signals

    -Methods:
        - read_channels()
        - get_raw_signal()
    """

    def __init__(self, filepath, plt_signals):
        self.filepath = filepath
        self.plt_signals = plt_signals
        self.logger = logging.getLogger('asammdf')
        self.logger.setLevel(50)

    def read_channels(self):
        """Extract the processable signal objects from the MDF file and
        make a narrowed asammdf.MDF object what contains only the signals
        what we need.

        -Parameters:
            - None

        -Returns:
            - new_mdf: [asammdf.MDF] an MDF object which contains only the \
                       data what is needed and successfully found in the \
                       input file
            - mdf_channels: [list] list of strings with the signal names of \
                            the new_mdf
        """
        # mdf_data is an asammdf.MDF object which contains all data from the
        # input file.

        try:
            mdf_data = asammdf.MDF(name=self.filepath)
            # signals is a list containing only the
            # neccessary data for the conversion (name, samples, timestamps)
            signals = []
            # the mdf_channels will containing the names of the signals at the
            # new_mdf
            mdf_channels = []
            # handler is an instance of SignalNameHandler class. It's needed
            # because of the signal names of the plt are not every time in the
            # same form as in the MDF file
            handler = SignalNameHandler(mdf_data)
            for channel in list(set(self.plt_signals)):
                signalmap = handler.find(channel)
                found = False
                signal_index = 0
                # If the handler finds more than one variation of the signal names
                # with insecure characters like space, underscore or colon we
                # should looking for all of the existing signals. Amongst each
                # signal we will use the very first with valid data.
                if signalmap:
                    names = list(signalmap.keys())
                    while (not found) and (len(names) > signal_index):
                        realname = names[signal_index]
                        locations = signalmap[realname]
                        location_index = 0
                        # The MDF standard allows us to store signal data with the
                        # same signal name in different signal groups.
                        # We will use the first valid data amongst them.
                        if locations:
                            while (not found) and (len(locations) > location_index):
                                try:
                                    groupindex, channelindex = locations[location_index]
                                    raw_signal = self.get_signal(mdf_data,
                                                                 realname,
                                                                 channel,
                                                                 groupindex,
                                                                 channelindex)
                                    if raw_signal:
                                        if (raw_signal.samples.dtype == numpy.float32
                                            or raw_signal.samples.dtype == numpy.float16):
                                            raw_signal.samples.astype(numpy.float64)
                                        signals.append([raw_signal.name, raw_signal.samples, raw_signal.timestamps])
                                        mdf_channels.append(channel)
                                        found = True
                                except Exception as signalerror:
                                    messages.messages.ErrorMessages("corrupted_signal", (channel,
                                                                                         groupindex,
                                                                                         channelindex,
                                                                                         signalerror))
                                location_index += 1
                        signal_index += 1
        except MemoryError:
            messages.messages.ErrorMessages("error_memory_limit")
            signals = None
            mdf_channels = None
        except Exception as error:
            messages.messages.ErrorMessages("asammdf_readerror", error)
            signals = None
            mdf_channels = None
        mdf_data.close()
        return signals, mdf_channels

    def get_signal(self, mdf_data, realname, channel, groupindex, channelindex):
        result = None
        requested_channel = [(realname, groupindex, channelindex)]
        raw_signal, processable = self.select_signal(mdf_data,
                                                     requested_channel,
                                                     channel,
                                                     raw=False,
                                                     ignore_value2text_conversions = False)
        if not processable:
            raw_signal, processable = self.select_signal(mdf_data,
                                                         requested_channel,
                                                         channel,
                                                         raw=False,
                                                         ignore_value2text_conversions = True)
            if not processable:
                raw_signal, processable = self.select_signal(mdf_data,
                                                             requested_channel,
                                                             channel,
                                                             raw=True)
                if not processable:
                    raw_signal = None
        if raw_signal:
            raw_signal.conversion = None
        if processable:
            result = raw_signal
        return result

    def select_signal(self, mdf_data, requested_channel, channel, raw, ignore_value2text_conversions = False):
        result = None
        processable = False
        channellist = mdf_data.select(requested_channel, raw=raw, ignore_value2text_conversions=ignore_value2text_conversions)
        for raw_signal in channellist:
            if (len(raw_signal.samples) > 0 and len(raw_signal.timestamps) > 0 and len(raw_signal.samples) == len(raw_signal.timestamps)):
                raw_signal.name = channel
                if self.is_processable(raw_signal.samples):
                    result = raw_signal
                    processable = True
        return result, processable

    def is_processable(self, samples):
        processable = True
        for value in samples:
            try:
                float(value)
                break
            except ValueError:
                processable = False
                break
        return processable


class SignalNameHandler:
    """The SignalNameHandler class provides a solution for insecure signal
    names (which conatains space, colon, underscore) in the input file or PLT

    -Init:
        - mdf_data: an asammdf.MDF object to search for the possible location\
                    of the signals

    -Methods:
        - find(pltname)
        - isinsecure(pltname)
        - insecurename(pltname)
    """
    def __init__(self, mdf_data):
        self.insecure_characters = [' ', ':']
        self.replacement_characters = ['_']
        self.basicpattern = r'(\s|\:|\_)'
        self.mdf_data = mdf_data
        self.channels = list(self.mdf_data.channels_db.keys())

    def find(self, pltname):
        """The main function of the SignalNameHandler class.
        In case of a secure name is equivalent with the asammdf.MDF.whereis()
        function.
        If the PLT signal is an insecure name and the asammdf.whereis() did not
        found in the measurement file we will looking for signal names in the
        mdf object with a regular expression which covers all variation of the
        signal name with insecure characters.

        -Parameters:
            - pltname: [string] the name of the signal which was defined in\
                       the plt

        -Returns:
            - signalmap: [dictionary] key(s): signal names in the plt\
                         (string),  value(s): [dictionary] key(s): signal\
                         names (string) in the mdf object,  value(s): tuple\
                         with two integer, the first is the groupindex the\
                         second is the channelindex
         """
        signalmap = {}
        realname = pltname
        locations = self.mdf_data.whereis(pltname)
        if locations == ():
            insecure = self.isinsecure(pltname)
            if not insecure:
                signalmap[realname] = locations
                return signalmap
            signalmap = self.insecurename(pltname)
            return signalmap
        signalmap[realname] = locations
        return signalmap

    def isinsecure(self, pltname):
        """A simple function what can decide if a signal name from the plt is
        insecure or not.

        -Parameters:
            - pltname: [string] the name of the signal which was defined in\
                       the plt

        -Returns:
            - [boolean]: True if the signal name is insecure False if the \
                         signal name is secure.
        """
        for character in self.replacement_characters:
            if character in pltname:
                return True
        return False

    def insecurename(self, pltname):
        """The function finds all possible signalname in case of insecure
        signalnames.

        -Parameters:
            - pltname: [string] the name of the signal which was defined in\
                       the plt

        -Returns:
            - signalmap: [dictionary] key(s): signal names (string) in the \
                         mdf object,  value(s): tuple with two integer, the \
                         first is the groupindex the second is the channelindex
        """
        signalmap = {}
        pltname_generalized = pltname.replace('\\', '\\\\').replace('_', self.basicpattern)
        pattern = re.compile(pltname_generalized)
        signalname_filter = filter(pattern.match, self.channels)
        signalnames = []
        for signalname in signalname_filter:
            signalnames.append(signalname)
        for signalname in signalnames:
            location = self.mdf_data.whereis(signalname)
            signalmap[signalname] = location
        return signalmap
