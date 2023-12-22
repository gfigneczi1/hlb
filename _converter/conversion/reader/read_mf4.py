"""
@author: fse2bp
"""

import logging
from _converter.libs.mf4reader import mf4reader


class ReadMF4:
    """The ReadMF4 class reads the MF4 file and contructs a signals Py list
    which contains the signals what is defined in the plt_signals.

    -Init:
        - filepath: the absolute path of an MF4 file
        - plt_signals: list of strings with the name of the needed signals

    -Methods:
        - read_channels()
    """

    def __init__(self, filepath, plt_signals):
        self.filepath = filepath
        self.plt_signals = plt_signals
        self.logger = logging.getLogger('mf4reader')
        self.logger.setLevel(50)

    def read_channels(self):
        """Extract the processable signal objects from the MF4 file and
        make signals Py list what contains only the signals
        what we need.

        -Parameters:
            - None

        -Returns:
            - signals: [list] an Python list which contains the \
                       signals what is needed and successfully found in the \
                       input file
            - mf4_channels: [list] list of strings with the signal names of \
                            the signals Py list
        """

        # mf4_data is an mf4reader object which contains the needed data from the
        # input file.
        if True: #try:
            mf4_data = mf4reader(filepath=self.filepath)
            # signals is a Py list containing only the
            # neccessary data for the conversion (name, samples, timestamps)
            signals = []
            # the mf4_channels will containing the names of the signals at the
            # new_mdf
            mf4_channels = []
            for channel in list(set(self.plt_signals)):
                channel_data = mf4_data.getSignalValue(channel)
                # signal = [Name, Value, Time]
                signal = [channel, channel_data['Value'], channel_data['Time']]
                signals.append(signal)
                mf4_channels.append(channel)

        #except MemoryError:
        #    messages.messages.ErrorMessages("error_memory_limit")
        #    signals = None
        #    mf4_channels = None
        #except Exception as error:
        #    messages.messages.ErrorMessages("mf4reader_readerror", error)
        #    signals = None
        #    mf4_channels = None

        mf4_data.close()
        return signals, mf4_channels
