/*
                BCMatch4Graphs
     Bootstrap AMG based on Compatible weighted Matching for Graphs, version 1.0
    (C) Copyright 2021
                       Pasqua D'Ambra         IAC-CNR, IT
                       Panayot S. Vassilevski Portland State University, OR USA

  Redistribution and use in source and binary forms, with or without
  modification, are permitted provided that the following conditions
  are met:
    1. Redistributions of source code must retain the above copyright
       notice, this list of conditions and the following disclaimer.
    2. Redistributions in binary form must reproduce the above copyright
       notice, this list of conditions, and the following disclaimer in the
       documentation and/or other materials provided with the distribution.
    3. The name of the BCMatch4Graphs group or the names of its contributors may
       not be used to endorse or promote products derived from this
       software without specific written permission.

  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
  ``AS IS'' AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
  TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
  PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE BCMATCH4GRAPHS GROUP OR ITS CONTRIBUTORS
  BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
  CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
  SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
  INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
  CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
  ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
  POSSIBILITY OF SUCH DAMAGE.
*/

#include "bcm_metadata_exp.h"

mge_Metadata* mge_MetadataCreate(int num_verts, int num_properties) {
    mge_Metadata* metadata = (mge_Metadata*) calloc(1, sizeof(mge_Metadata));
    metadata -> num_verts = num_verts;
    metadata -> num_properties = num_properties;
    metadata -> metadata_matrix = (int*) calloc(num_verts * num_properties, sizeof(int));
    return metadata;
}

int mge_MetadataGetLabel(mge_Metadata* metadata, int vertex, int property)
{
    if (vertex > (metadata -> num_verts) || property > (metadata -> num_properties)) {
        printf("Warning! Metadata retrieval out of bounds!\n:");
    }

    int index = (metadata -> num_properties) * vertex + property;
    return (metadata -> metadata_matrix)[index];
}

int mge_MetadataDestroy(mge_Metadata* metadata) {
  int  ierr=0;
  if (metadata)
    {
      free(metadata -> metadata_matrix);
      free(metadata);
    }
}
