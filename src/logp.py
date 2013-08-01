# -*- coding: utf-8 -*-

##    Description    eTOXlab model template
##                   
##    Authors:       Manuel Pastor (manuel.pastor@upf.edu) 
##
##    Copyright 2013 Manuel Pastor
##
##    This file is part of eTOXlab.
##
##    eTOXlab is free software: you can redistribute it and/or modify
##    it under the terms of the GNU General Public License as published by
##    the Free Software Foundation version 3.
##
##    eTOXlab is distributed in the hope that it will be useful,
##    but WITHOUT ANY WARRANTY; without even the implied warranty of
##    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##    GNU General Public License for more details.
##
##    You should have received a copy of the GNU General Public License
##    along with eTOXlab.  If not, see <http://www.gnu.org/licenses/>.

from rdkit.Chem import Descriptors
from rdkit import Chem


def computeLogP (mol):

    lp = []
    try:
        suppl = Chem.SDMolSupplier(mol)
        mi = suppl.next()

        if mi is None:
            return (False, 'wrong input format')

        lp.append(Descriptors.MolLogP(mi))
        
    except:
        return (False, 'wrong input format')

    if len(lp) == 0:
        return (False,'error in logP computation')
    else:
        return (True, lp)  
