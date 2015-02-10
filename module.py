#!/opt/local/bin/ipython
# -*- coding: utf-8 -*-
from logging_init import *
import __main__ as main
main_file = main.__file__

class Personne(object):
    def __init__(self, nom):
        # Logging
        self.loggerm = init_log(__name__)
        self.nom = nom
        self.loggerm.info(nom)




    def presente(self):
        self.loggerm.info(self.nom)
        self.loggerm.warn('Création!!!!')



def myFunction():
    print('ma fonction')
