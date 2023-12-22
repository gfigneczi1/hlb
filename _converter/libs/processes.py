"""
Created on Tue Jan 22 11:45:03 2019
@author: jez2bp
"""

import getpass
import time
import logging
from requests.exceptions import ChunkedEncodingError
from urllib3.exceptions import ProtocolError
from appdata.buildinfo import Buildinfo
from conversion.convert_main import MainConversion
from clientapi.sessions import Sessions
from clientapi.downloads import downloads
from clientapi.dataqueue import Dataqueue
from clientapi.processes import Processes

class ConvertProcesses:
    """The class provides functions for the running of the conversion tool.
    Handles the command line arguments and the passing between the functions.
    """
    def __init__(self, config, settings, vehicle_profile, threadswitch, hookprogress=None, hookcredentials=None):
        logging.getLogger("requests").setLevel(50)
        logging.getLogger("urllib3").setLevel(50)
        self.settings = settings
        self.config = config
        self.vehicle_profile = vehicle_profile
        self.threadswitch = threadswitch
        self.hookprogress = hookprogress
        if hookcredentials is not None:
            self.hookcredentials = hookcredentials
        else:
            self.hookcredentials = self.consoleCredentials
        build = Buildinfo()
        self.build = build.read()
        self.onlinesession = None
        self.username = None
        # Options
        if self.settings["evals"] or self.settings["uploadstate"] or self.settings["loginstate"]:
            try:
                self.onlinesession = self.createSession()
            except (ConnectionResetError, ChunkedEncodingError, ProtocolError):
                logging.error('Connection is not possible!')
                self.onlinesession = None
                self.settings["evals"] = []
                self.settings["uploadstate"] = False
                self.settings["loginstate"] = False
        self.handling_contradictory_parameters()

    def run(self):
        """The method starts the conversion and the uploading"""
        self.settings["evalsignals"] = {}
        if self.settings["evals"] and self.onlinesession.sessionid is not None:
            self.settings["evalsignals"] = self.listProcesssignals(self.settings["evals"])
        converter = MainConversion(self.settings, self.vehicle_profile, self.threadswitch, self.hookprogress)
        files_for_upload = converter.start()
        if self.settings["uploadstate"] and files_for_upload is not None:
            self.uploadMeasure(files_for_upload)
        if self.onlinesession is not None:
            self.onlinesession.logout()
        logging.info("Completed")

    def uploadMeasure(self, files_for_upload):
        """Uploads the converted files into the Online KPI Tool"""
        dataqueue = Dataqueue(config=None, apiurl=self.config.get("network", "apiurl"))
        dataqueue.sessionid = self.onlinesession.sessionid
        logging.info("Uploading...")
        dataqueue.pushData(files_for_upload)
        if int(dataqueue.lastanswer.status_code) != 200:
            logging.error("Something went wrong while uploading")
            logging.error("Status code: " + str(dataqueue.lastanswer.status_code))

    def createSession(self):
        """Creates a session for the cloud connection needed functions"""
        result = Sessions(config=None, apiurl=self.config.get("network", "apiurl"))
        userid = self.tokenLogin(result)
        if userid is None:
            userid = self.login(result)
        if userid is not None:
            self.checkVersion(result, self.build)
        return result

    def checkVersion(self, session, build):
        download = downloads(config=None, apiurl=self.config.get("network", "apiurl"))
        download.sessionid = session.sessionid
        freshbuild = download.checkversion('win-conv4kpi')
        if int(freshbuild) > build:
            logging.warning("Your current Conv4KPI build is outdated! Please use the newest version! visit http://kpi.apps.de1.bosch-iot-cloud.com")

    def consoleCredentials(self, username):
        username = str(input("Username: "))
        password = getpass.getpass("Password: ")
        return len(password)>0, username, password
        
    def tokenLogin(self, session):
        """Login the session with authentication token if its available"""
        result = None
        token = self.config.get("network", "token")
        if token:
            #userid
            result = session.tokenLogin(token)
        return result

    def login(self, session):
        """Login the session with username and password if the token expired or
        not available"""
        result = None
        logging.info("You are using at least one option what requires cloud connection!")
        logging.info("Please login into the Online KPI Tool account or register...")
        logging.info("http://kpi.apps.de1.bosch-iot-cloud.com")
        resp = False
        while not resp:
            state, self.username, password = self.hookcredentials(self.username)
            if state:
                #userid
                result = session.login(self.username, password)
                if result is None:
                    logging.error("Authentication failed! Try again!")
                    time.sleep(3)
                else:
                    token = session.requireToken()
                    self.config.set('network', 'token', token)
                    self.config.save()
                    resp = True
            else:
                #if state is false, the user has canceled the input, so dont repeat
                resp = True
        return result

    def listProcesssignals(self, evaluation):
        """Gets the necessary signal names for the given --eval parameter"""
        result = {}
        logging.info("Requesting signal names for signal evaluation checking...")
        processes = Processes(config=None, apiurl=self.config.get("network", "apiurl"))
        processes.sessionid = self.onlinesession.sessionid
        for eval_type in evaluation:
            result[eval_type] = processes.listProcesssignals(eval_type)
        logging.info("Completed")
        return result

    def handling_contradictory_parameters(self):
        """Handles the contradictory parameters and if its needed change the option argument value"""
        if self.onlinesession is None and self.settings["uploadstate"]:
            logging.error("Upload is not possible because authentication failed!")
            self.settings["uploadstate"] = False
        if self.onlinesession is None and self.settings["loginstate"]:
            logging.error("Login is not possible because authentication failed!")
            self.settings["loginstate"] = False
        if self.onlinesession is None and self.settings["evals"]:
            logging.error("Signal evaluation checking is not possible because authentication failed!")
            self.settings["evals"] = None
