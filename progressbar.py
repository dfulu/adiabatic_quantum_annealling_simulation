# -*- coding: utf-8 -*-
import sys
import time

# Print iterations progress
def printProgress (iteration, total, prefix = '', suffix = '', decimals = 1, barLength = 100, PROGRESS_IT_NUMBER = None):
    """
    Call in a loop to create terminal progress bar
    @params:
        iteration   - Required  : current iteration (Int)
        total       - Required  : total iterations (Int)
        prefix      - Optional  : prefix string (Str)
        suffix      - Optional  : suffix string (Str)
        decimals    - Optional  : positive number of decimals in percent complete (Int)
        barLength   - Optional  : character length of bar (Int)
    """
    if iteration == 0: 
        global t_i
        t_i= time.time()
    else:
        if suffix== 'est_time':
            suffix = '|Est. Time Remaining : %.1f'%((time.time()-t_i)*(float(total)/iteration-1)) + ' '*10
    
    formatStr       = "{0:." + str(decimals) + "f}"
    percents        = formatStr.format(100 * (iteration / float(total)))
    filledLength    = int(round(barLength * iteration / float(total)))
    bar             = '█' * filledLength + '-' * (barLength - filledLength)
    sys.stdout.write('\r%s |%s| %s%s %s' % (prefix, bar, percents, '%', suffix)),
    sys.stdout.flush()
    if iteration == total: 
        dt = time.time()-t_i
        rt =  'Runtime: %sh %sm %ss'%(int(dt)/3600,(int(dt)/60)%60, '{0:.2f}'.format(dt%60))
        if PROGRESS_IT_NUMBER != None:
            rt += ' - Completed run: %s'%PROGRESS_IT_NUMBER
        print'\r'+rt+' '*(-len(rt)+len('\r%s |%s| %s%s %s' % (prefix, bar, percents, '%', suffix)))
    
def printProgress2part (finishedfirstpart, iteration, total, prefix = '', suffix = '', decimals = 1, barLength = 100, PROGRESS_IT_NUMBER = None):
    """
    Call in a loop to create terminal progress bar
    @params:
        iteration   - Required  : current iteration (Int)
        total       - Required  : total iterations (Int)
        prefix      - Optional  : prefix string (Str)
        suffix      - Optional  : suffix string (Str)
        decimals    - Optional  : positive number of decimals in percent complete (Int)
        barLength   - Optional  : character length of bar (Int)
    """
    global t_i_firstpart
    global t_i
    if not finishedfirstpart:
        if iteration == -1: 
            
            t_i_firstpart= time.time()
        else:
            if suffix== 'est_time':
                suffix = '|Est. Time Remaining : %.1f'%((time.time()-t_i_firstpart)*(float(total)/iteration-1)) + ' '*10
        
        
        formatStr       = "{0:." + str(decimals) + "f}"
        percents        = formatStr.format(100 * (iteration / float(total)))
        filledLength    = int(round(barLength * iteration / float(total)))
        bar             = '-'*10+'|'+'+' * filledLength + '-' * (barLength - filledLength)
        sys.stdout.write('\r%s |%s| %s%s %s' % (prefix, bar, percents, '%', suffix)),
        sys.stdout.flush()
        if iteration == total: 
            dt = time.time()-t_i
            rt =  'Runtime: %sh %sm %ss'%(int(dt)/3600,(int(dt)/60)%60, '{0:.2f}'.format(dt%60))
            if PROGRESS_IT_NUMBER != None:
                rt += ' - Completed run: %s'%PROGRESS_IT_NUMBER
            print'\r'+rt+' '*(-len(rt)+len('\r%s |%s| %s%s %s' % (prefix, bar, percents, '%', suffix)))


    
    elif finishedfirstpart:
        if iteration == 0: 
            
            t_i= time.time()
        else:
            if suffix== 'est_time':
                suffix = '|Est. Time Remaining : %.1f'%((time.time()-t_i)*(float(total)/iteration-1)) + ' '*10
        
        formatStr       = "{0:." + str(decimals) + "f}"
        percents        = formatStr.format(100 * (iteration / float(total)))
        filledLength    = int(round(barLength * iteration / float(total)))
        bar             = '#'*10+'|'+'+' * filledLength + '-' * (barLength - filledLength)
        sys.stdout.write('\r%s |%s| %s%s %s' % (prefix, bar, percents, '%', suffix)),
        sys.stdout.flush()
        if iteration == total: 
            dt = time.time()-t_i
            rt =  'Runtime: %sh %sm %ss'%(int(dt)/3600,(int(dt)/60)%60, '{0:.2f}'.format(dt%60))
            if PROGRESS_IT_NUMBER != None:
                rt += ' - Completed run: %s'%PROGRESS_IT_NUMBER
            print'\r'+rt+' '*(-len(rt)+len('\r%s |%s| %s%s %s' % (prefix, bar, percents, '%', suffix)))


# 
# Sample Usage
# 

if __name__ == "__main__":
    from time import sleep
    # make a list
    items = list(range(0, 100))
    l     = len(items)-1

    # Initial call to print 0% progress
    for i in range(4):
        PROGRESS_IT_NUMBER = i+1
        printProgress2part(False,-1, l, prefix = 'Progress:', suffix = 'est_time', barLength = 20, PROGRESS_IT_NUMBER = i+1)
        sleep(2)
        for item in items:
            # Do stuff...
            sleep(0.01)
            # Update Progress Bar
            printProgress2part(True,item, l, prefix = 'Progress:', suffix = 'est_time', barLength = 20, PROGRESS_IT_NUMBER = i+1)