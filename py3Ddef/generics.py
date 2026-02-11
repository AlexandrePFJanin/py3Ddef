# -*- coding: utf-8 -*-
"""
@author: Alexandre JANIN
@aim:    Generic routines
"""

# External dependencies:
from termcolor import colored


# ----------------- FUNCTIONS -----------------


def im(textMessage,pName,verbose,error=False,warn=False,structure=False,end=False):
    """Print verbose internal message. This function depends on the
    argument verbose. If verbose, then the message will be displayed
    in the terminal.
    
    Args:
        textMessage = str, message to display
        pName = str, name of the subprogram
        verbose = bool, condition for the verbose output
    """
    if verbose and not error:
        if not warn:
            msgc = None
        else:
            msgc = 'yellow'
        if structure:
            if end:
                print(colored('>> '+pName+'| ','blue')+colored('--- ','magenta')+colored(textMessage, msgc))
            else:
                print(colored('>> '+pName+'| ','blue')+colored(' : ','magenta')+colored(textMessage, msgc))
        else:
            print(colored('>> '+pName+'| ','blue')+colored(textMessage, msgc))
    if error:
        #print error message
        print(colored('>> '+pName+'| --- ----- ---','red'))
        print(colored('>> '+pName+'| --- ERROR ---','red'))
        print(colored('>> '+pName+'| --- ----- ---','red'))
        print(colored('>> '+pName+'| '+textMessage,'red'))
        raise AssertionError()