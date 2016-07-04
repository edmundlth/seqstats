#Overview
---
Seqstats is a tool intended to compute and report various statistics of a DNA sequence (embeded in a reference sequence)and its "context". A context is used here to mean the flanking sequences of a fixed predefined size before and after a sequence with respect to a reference sequence. 

Current supported statistics:
    1. Sequence entropy
    2. Context sequences entropy

The goal is to provide quick statistical information to the users. As a suggestion for future improvement, statistical significance test of these statistics against random sampling from the reference sequence can be implemented.


#Licence
---
This program is released as open source software under the terms of [MIT License](https://github.com/edmundlth/seqstats/blob/master/LICENCE)

#Installing
---

Seqstats can be installed using pip in a variety of ways (% indicates the command line prompt):

1. Inside a [virtual environment](https://virtualenv.pypa.io/en/stable/):  
`% virtualenv seqstats_env`  
`% source seqstats_env/bin/activate`  
`% pip install -U /path/to/seqstats-py`  
2. Into the global package database for all users:  
`% pip install -U /path/to/seqstats-py`
3. Into the user package database (for the current user only):  
`% pip install -U --user /path/to/seqstats-py`

#General Behaviour
---

#Usage
---
##Help message



#Error Handling
---
##Exit Status Values


#Testing
---

#Bugs
---
Filed at [Issue Tracker](https://github.com/edmundlth/seqstats/issues)