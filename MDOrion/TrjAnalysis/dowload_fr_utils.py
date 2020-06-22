# (C) 2020 OpenEye Scientific Software Inc. All rights reserved.
#
# TERMS FOR USE OF SAMPLE CODE The software below ("Sample Code") is
# provided to current licensees or subscribers of OpenEye products or
# SaaS offerings (each a "Customer").
# Customer is hereby permitted to use, copy, and modify the Sample Code,
# subject to these terms. OpenEye claims no rights to Customer's
# modifications. Modification of Sample Code is at Customer's sole and
# exclusive risk. Sample Code may require Customer to have a then
# current license or subscription to the applicable OpenEye offering.
# THE SAMPLE CODE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
# EXPRESS OR IMPLIED.  OPENEYE DISCLAIMS ALL WARRANTIES, INCLUDING, BUT
# NOT LIMITED TO, WARRANTIES OF MERCHANTABILITY, FITNESS FOR A
# PARTICULAR PURPOSE AND NONINFRINGEMENT. In no event shall OpenEye be
# liable for any damages or liability in connection with the Sample Code
# or its use.

import fileinput

import re

import tempfile

import os

import zipfile

from orionclient.types import ShardCollection

resource_url = "/api/v1/storage/collections/{}/shards/{}/download/"
page_url = "/app/main/jobreport2/{}/{}/{}"


def replace_links(filename, link_map):
    def resource_to_file(match):
        shard_id = int(match.group(2))
        return link_map[shard_id]

    def page_to_file(match):
        shard_id = int(match.group(3))
        return link_map[shard_id]

    match_number = "([0-9]+)"
    shard_regex = resource_url.format(match_number, match_number)
    page_regex = page_url.format(match_number, match_number, match_number)
    with fileinput.FileInput(filename, inplace=True, backup='.bak') as file:
        for line in file:
            print(re.sub(page_regex, page_to_file,
                         re.sub(shard_regex, resource_to_file, line)))


def download_report(collection_id, session, zip_filename="floe_report.zip"):

    collection = session.get_resource(ShardCollection, collection_id)
    index_pages = collection.metadata["index_pages"]
    index_ids = list(page["id"] for page in index_pages)
    shards = list(collection.list_shards())
    index_counter = 1
    other_page_counter = 1
    shard_file_map = {}

    filenames = []
    with tempfile.TemporaryDirectory() as temp_dir:
        os.mkdir(os.path.join(temp_dir, "resources"))
        for shard in shards:
            if shard.id in index_ids:
                name = "index_" + str(index_counter) + ".html"
                index_counter = index_counter + 1
            else:
                name = os.path.join("resources", "resource_" + str(other_page_counter) + ".html")
                other_page_counter = other_page_counter + 1

            fullpath = os.path.join(temp_dir, name)
            shard.download_to_file(fullpath)
            shard_file_map[shard.id] = name
            filenames.append([fullpath, name])

        for fullpath, fname in filenames:
            replace_links(fullpath, shard_file_map)

        print("Creating Floe Report zip file...")
        zf = zipfile.ZipFile(zip_filename, mode='w')
        try:
            for fullpath, fname in filenames:
                zf.write(fullpath, fname)
        finally:
            zf.close()

        return zip_filename
