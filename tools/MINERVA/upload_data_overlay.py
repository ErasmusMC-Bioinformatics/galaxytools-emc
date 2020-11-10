#!/usr/bin/env python
"""
MINERVA API version: 14.0
Script purpose: An example of MINERVA API calls uploading overlays
Date: 23/03/2020
Author: Piort Gawron (piotr.gawron@uni.lu)
MANUAL: https://minerva.pages.uni.lu/doc/api/14.0/projects/
GALAXY: Willem de Koning 10/11/2020
##################################################
"""

import requests
import os
import json
import argparse


def upload_overlay(login, password, overlay_file, api_url, project_id):
    session = requests.Session()

    # LOGIN to PD map instance with your credentials
    login_request = session.post(api_url + "/doLogin",
                                 data={'login': login, 'password': password})
    print(login_request.text)
    print(session.cookies.get_dict())

    # UPLOAD OVERLAY FILE TO THE SYSYETM
    stat_info = os.stat(overlay_file)
    with open(overlay_file) as f:
        file_content = f.read()

    # Before adding an overlay in MINERVA, data overlay file must be uploaded
    # to the instance
    # (MANUAL: https://minerva.pages.uni.lu/doc/api/14.0/files/)
    # Allocate memory in the system for 'overlay_file', which length is
    # 'stat_info.st_size' bytes
    create_file_request = session.post(api_url+"/files/",
                                       data={'filename': overlay_file,
                                             'length': stat_info.st_size})

    # Get the information about the uploaded file:
    # the 'id' is necessary to upload the file's content
    content = json.loads(create_file_request.text)

    file_id = content["id"]
    # Upload file's content to the instance
    upload_content_request = session.post(api_url + "/files/" +
                                          str(file_id) + ":uploadContent",
                                          data=file_content)

    # CREATE DATA OVERLAY FROM UPLOADED FILE
    create_overlay_request = session.post(api_url + "/projects/" +
                                          project_id + "/overlays/",
                                          data={"fileId": file_id,
                                                "filename": overlay_file,
                                                "name":
                                                "auto-overlay-simple" +
                                                str(file_id),
                                                "description":
                                                "there is a data overlay",
                                                "googleLicenseConsent":
                                                "true"})

    print("Overlay Added: ", create_overlay_request.text)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='An example of MINERVA API' +
                                                 'calls uploading overlays')
    parser.add_argument('--login', help='Login name', required=True)
    parser.add_argument('--password', help='Password', required=True)
    parser.add_argument('--overlay', help='Overlay file', required=True)
    parser.add_argument('--api_url', help='API url', required=True)
    parser.add_argument('--project_id', help='Project ID', required=True)
    args = parser.parse_args()

    login = args.login
    password = args.password
    overlay_file = args.overlay
    api_url = args.api_url
    project_id = args.project_id

    upload_overlay(login, password, overlay_file, api_url, project_id)

