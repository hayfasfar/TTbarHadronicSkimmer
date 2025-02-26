# ttbaranalysis.py

from coffea import util
from coffea.nanoevents import NanoAODSchema, BaseSchema
import coffea.processor as processor

import itertools
import argparse
import time
from datetime import date
import json
import os

from dask.distributed import Client, performance_report

import warnings
warnings.filterwarnings("ignore")

default_datastets = ['JetHT', 'TTbar', 'QCD']
default_signals = ['RSGluon', 'ZPrime10', 'ZPrime30', 'ZPrimeDM', 'ZPrime1']

from ttbarprocessor import TTbarResProcessor
from python.functions import printTime, makeSaveDirectories






if __name__ == "__main__":
    
    tic = time.time()

    

    # directory where output files are saved
    #savedir = f'outputs/ttagSF/' 
    savedir = f'outputs/ttbar16/' 
    
    parser = argparse.ArgumentParser(
                    prog='ttbaranalysis.py',
                    description='Run ttbarprocessor',
                    epilog='for a test of MC QCD, MC TTbar, JetHT run  "python ttbaranalysis.py --test"')
    
    # datasets to run
    parser.add_argument('-d', '--dataset',
                        choices=['JetHT', 'QCD', 'TTbar', 'ZPrime1', 'ZPrime10', 'ZPrime30', 'ZPrimeDM', 'RSGluon'], 
                        default=default_datastets,
                        action='append'
                       )
    
    parser.add_argument('--iov', choices=['2016APV', '2016', '2017', '2018'], default='2016')
    parser.add_argument('--signals', action='store_true', help='run only signal samples')
    
    
    
    
    # choose specific eras, pt bins, mass points
    parser.add_argument('--era', choices=['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H'], action='append', default=[], help='--era A --era B --era C for multiple eras, runs all eras if not specificed')
    parser.add_argument('-p', '--pt', choices=['700to1000', '1000toInf'], action='append', default=[], help='pt bins for TTbar datasets')
    parser.add_argument('-m', '--mass', action='append', default=[], help='mass points for signal')

    # analysis options
    parser.add_argument('--blind', action='store_true', help='process 1/10th of the data')
    parser.add_argument('--bkgest', choices=['2dalphabet', 'mistag'], default=None)
    parser.add_argument('--toptagger', choices=['deepak8', 'cmsv2'], default='deepak8')
    parser.add_argument('-r', '--redirector', default='root://cmsxrootd.fnal.gov/')
    parser.add_argument('--ttagWP', choices=['loose', 'medium', 'tight'], default='medium')
    parser.add_argument('--btagger', choices=['deepcsv', 'csvv2'], default='deepcsv')
    parser.add_argument('--ht', choices=['1400', '950'], default='1400')
    parser.add_argument('--noSyst', action='store_true', help='run without systematics')

    # run options
    parser.add_argument('--dask', action='store_true')
    parser.add_argument('--env', choices=['casa', 'lpc', 'winterfell', 'C', 'L', 'W'], default='lpc')
    parser.add_argument('--test', action='store_true')
    parser.add_argument('-n', '--nocluster', action='store_true', help='use client=Client() if LPCCondorCluster is slow')

    args = parser.parse_args()
    
    # remove defaults if --dataset given
    if len(args.dataset) > len(default_datastets): 
        args.dataset = args.dataset[len(default_datastets):]
        
    if args.signals:
        args.dataset = default_signals
    
    if args.dask and (args.env == 'lpc' or args.env == 'L'):
        from lpcjobqueue import LPCCondorCluster
    
    
    ##### parameters #####

    samples = args.dataset
    IOV = args.iov
    useDeepAK8 = True if (args.toptagger == 'deepak8') else False
    useDeepCSV = True if (args.btagger == 'deepcsv') else False
    htCut = 1400.0 if (args.ht == '1400') else 950.0
    dask_memory = '3GB' # priority decreases for >2GB memory
    chunksize_dask = 100000
    chunksize_futures = 10000
    maxchunks = 10 if args.test else None
    
    
    
    # transfer function parameters
    #params = json.load(open(f'data/corrections/params_{IOV}.json'))


        
    
    ##### systematics and analysis categories #####
    
    ############################################################################################################
    #                                                                                                          #   
    # systematics options: nominal, jes, jer, pileup, pdf, q2, btag, prefiring                                 # 
    #                                                                                                          #   
    # top tag category options:                                                                                #
    #                                                                                                          #
    # at        antitag - jet0 is antitagged                                                                   #            
    # pret      pretag  - jet0 is top tagged, jet1 top tag is unknown                                          #                 
    # 0t        jet0 and jet1 are not top tagged                                                               #                
    # 1t        jet0 is top tagged and jet1 is not top tagged                                                  #
    # 2t        jet0 and jet1 are top tagged                                                                   #            
    # AT&P      antitag and probe - jet0 is antitagged and jet1 is top tagged                                  #               
    # >=1t      at least 1 top tag in event                                                                    #           
    # >=2t      at least 2 top tags in event                                                                   #            
    #                                                                                                          #             
    ############################################################################################################
    
    systematics = [
        'nominal',
        'jes',
        'jer',
        'pileup',
        'pdf',
        'q2',
        'ttag_pt1',
        'ttag_pt2',
        'ttag_pt3'
    ]
    
    if ('2016' in IOV) or ('2017' in IOV): systematics.append('prefiring')

    if args.bkgest == '2dalphabet': systematics.append('transferFunction')

     
    # make analysis categories 
    ttagcats = ["at", "2t"]
    ycats = ['cen', 'fwd']
    
    anacats = [ t+y for t,y in itertools.product( ttagcats, ycats) ]
    label_map = {i: label for i, label in enumerate(anacats)}
        
    
    
    # save analysis info #
    with open('out.log', 'w') as f:  
        print('\n'+date.today().isoformat(), file=f)
        print('\n------args------', file=f)
        for argname, value in vars(args).items(): print(argname, '=', value, file=f)
        print('----------------\n', file=f)
        print('categories =', label_map, file=f)
        print('\n', file=f)
        if not args.noSyst: print('systematics =', systematics, file=f)
     
    # display analysis info
    print('\n------args------')
    for argname, value in vars(args).items(): print(argname, '=', value)
    if not args.noSyst: print('systematics =', systematics)
    print('----------------\n')

    
    
    ##### root files #####
    if args.env == 'casa' or args.env == 'C': redirector = 'root://xcache/'
    elif args.env == 'winterfell' or args.env == 'W': redirector = '/mnt/data/cms/'
    else: redirector = redirector = args.redirector #'root://cmsxrootd.fnal.gov/' # default LPC
# #     redirector = 'root://xrootd-local.unl.edu:1094/'
#     # redirector = 'root://cmseos.fnal.gov:1094/'
#     redirector = 'root://cmseos.fnal.gov/'

    

    # if ('2018' in IOV):
        # redirector = 'root://eoscms.cern.ch//eos/cms/'
        # redirector = 'root://xrootd-local.unl.edu:1094/'
        # redirector = 'root://cmsxrootd.hep.wisc.edu:1094/' # Zprime1 2018
    # elif ('2017' in IOV):
        # redirector = 'root://cmsxrootd.hep.wisc.edu:1094/'
        # redirector = 'root://cms03.lcg.cscs.ch:1094/'
        # redirector = 'root://xrootd-local.unl.edu:1094/'
    # elif ('2016' in IOV):
    #     redirector = 'root://cmseos.fnal.gov/'




    
 
    jsonfiles = {
        "JetHT": 'data/nanoAOD/JetHT.json',
        "QCD": 'data/nanoAOD/QCD.json',
        "TTbar": 'data/nanoAOD/TTbar.json',
        "ZPrime1": 'data/nanoAOD/ZPrime1.json',
        "ZPrime10": 'data/nanoAOD/ZPrime10.json',
        "ZPrime30": 'data/nanoAOD/ZPrime30.json',
        "ZPrimeDM": 'data/nanoAOD/ZPrimeDM.json',
        "RSGluon": 'data/nanoAOD/RSGluon.json',
    }
    
    # directories and files for dask
    upload_to_dask = ['data', 'python', 'ttbarprocessor.py']



    # make output directories, save copy of processor, log changes and jobs run 
    if not os.path.exists(savedir):
        os.makedirs(savedir)
        os.makedirs(savedir+'logs/')
        os.makedirs(savedir+'scale/')
        os.makedirs(savedir+'twodalphabet/')
        os.popen('cp ttbarprocessor.py '+savedir+'logs/ttbarprocessor_'+date.today().isoformat().replace('-','')+'.py') 
        os.popen('cat out.log >> '+savedir+'logs/ttbarprocessor_diff.txt')
    
    else:
        for f in os.listdir(savedir+'logs/'):
            if 'ttbarprocessor' in f and 'py' in f:
                os.popen('cat out.log >> '+savedir+'logs/ttbarprocessor_diff.txt')
                print('diff ttbarprocessor.py '+savedir+'logs/'+f+' >> '+savedir+'logs/ttbarprocessor_diff.txt')
                os.popen('diff ttbarprocessor.py '+savedir+'logs/'+f+' >> '+savedir+'logs/ttbarprocessor_diff.txt')
                                
            
        if not os.path.exists(savedir):
            os.makedirs(savedir+'scale/')
        if not os.path.exists(savedir+'twodalphabet/'):
            os.makedirs(savedir+'twodalphabet/')


    # make image directories
    makeSaveDirectories(coffea_dir=savedir)
    

    
    
    ##### get fileset and run processor #####

    for sample in samples:
        
        # skipbadfiles = False if (('JetHT' in sample) or ('RSGluon' in sample) or ('ZPrime' in sample)) else True
        skipbadfiles = False
        #skipbadfiles = True

        # if 'ZPrime' in sample and 'ZPrime1' != sample:
        #         redirector = 'root://cmseos.fnal.gov/'
        # if 'ZPrime1' in sample:
        #         redirector = 'root://cmsxrootd.fnal.gov/'
        #     redirector = 'root://cmsxrootd.hep.wisc.edu:1094/'
        
        inputfile = jsonfiles[sample]
        files = []
        with open(inputfile) as json_file:
            
            
            subsections = args.era + args.mass + args.pt

            data = json.load(json_file)

            # select files to run over
            filedict = {}
            
            # if IOV split into subsections (eras, masses, pt bins)
            try: 
                data[IOV].keys()
                
                # if eras, pt bins, or mass points specified, add files individually
                if len(subsections) > 0:
                    for s in subsections:
                        if s in data[IOV].keys():
                            filedict[s] = data[IOV][s]
                        else:
                            print(f'{s} not in {sample} {IOV}')
                            
                            
                # if eras not specified, get all files for dataset
                else:
                    filedict = data[IOV]
                
            # no subsections in IOV    
            except:
                filedict[''] = data[IOV]
                

            # run uproot job
            for subsection, files in filedict.items():
                

                # add redirector; select file for testing
                files = [redirector + f for f in files]
                if args.test: files = [files[int(len(files)/2)]]
                    
                
                fileset = {sample: files}  
                                                 
                print(files[0])                  

                # coffea output file name
                subString = subsection.replace('700to', '_700to').replace('1000to','_1000to')
                if args.bkgest: subString += '_bkgest'
                    
                    

                if (args.toptagger == 'cmsv2') and (args.btagger == 'csvv2'): savedir = 'outputs/oldanalysis/'



                



                savefilename = f'{savedir}{sample}_{IOV}{subString}.coffea'
                if 'RSGluon' in sample:
                    subString = subString.replace(subsection, '')
                    savefilename = f'{savedir}{sample}{subsection}_{IOV}{subString}.coffea'
                elif 'ZPrime' in sample:
                    subString = subString.replace(subsection, '')
                    savefilename = f'{savedir}ZPrime{subsection}_{sample.replace("ZPrime","")}_{IOV}{subString}.coffea'
                print(f'running {IOV} {sample} {subsection}')

                if args.toptagger == 'cmsv2':
                    savefilename = savefilename.replace('.coffea', '_cmsv2.coffea')
                if args.btagger == 'csvv2':
                    savefilename = savefilename.replace('.coffea', '_csvv2.coffea')
                if args.ht == '950':
                    savefilename = savefilename.replace('.coffea', '_ht950.coffea')

                
                if args.blind: savefilename = savefilename.replace('.coffea', '_blind.coffea')
                if args.noSyst: savefilename = savefilename.replace('.coffea', '_noSyst.coffea')
                if args.test: savefilename = savefilename.replace('.coffea', '_test.coffea')
                
                
                
                # run using futures executor
                if not args.dask:

                    output, metrics = processor.run_uproot_job(
                        fileset,
                        treename="Events",
                        processor_instance=TTbarResProcessor(
                                                             iov=IOV,
                                                             bkgEst=args.bkgest,
                                                             noSyst=args.noSyst,
                                                             deepAK8Cut=args.ttagWP,
                                                             useDeepAK8=useDeepAK8,
                            
                                                             useDeepCSV=useDeepCSV,
                                                             htCut=htCut,
                                                             anacats=anacats,
                                                             systematics=systematics,
                                                             blinding=args.blind,
                                                             #rpf_params = params,

                                                            ),
                        executor=processor.futures_executor,
                        executor_args={
                                "skipbadfiles": skipbadfiles,
                                "savemetrics": True,
                                "schema": NanoAODSchema,
                                "workers":4
                                },
                        chunksize=chunksize_futures,
                        maxchunks=maxchunks,
                    )


                # run using dask
                else:
                    
                    if args.dask and (args.env == 'lpc' or args.env == 'L'):
                        
                        if args.nocluster:
                            cluster = None
                        else:
                            cluster = LPCCondorCluster(memory=dask_memory, transfer_input_files=upload_to_dask)
                            cluster.adapt(minimum=1, maximum=100)
                    else:
                        cluster = None
        

                    # files and directories for dask
                    upload_to_dask = [
                        'data',
                        'python',
                        'ttbarprocessor.py',
                    ]

                    
                    with Client(cluster) as client:

                        run_instance = processor.Runner(
                            metadata_cache={},
                            executor=processor.DaskExecutor(client=client, retries=12,),
                            schema=NanoAODSchema,
                            savemetrics=True,
                            skipbadfiles=skipbadfiles,
                            chunksize=chunksize_dask,
                            maxchunks=maxchunks,
                        )


                        if args.nocluster:
                            worker_toc = time.time()
                            print("Waiting for 4 workers...")
                            client.wait_for_workers(4)
                            worker_tic = time.time()
                            
                        else:
                            worker_toc = time.time()
                            print("Waiting for at least one worker...")
                            client.wait_for_workers(1)
                            worker_tic = time.time()
                        
                        print(f'time to wait for worker = {int(worker_tic - worker_toc)}s')

                        output, metrics = run_instance(fileset,
                                                      treename="Events",
                                                      processor_instance=TTbarResProcessor(
                                                          iov=IOV,
                                                          bkgEst=args.bkgest,
                                                          noSyst=args.noSyst,
                                                          deepAK8Cut=args.ttagWP,
                                                          useDeepAK8=useDeepAK8,
                                                          useDeepCSV=useDeepCSV,
                                                          htCut=htCut,
                                                          anacats=anacats,
                                                          systematics=systematics,
                                                          blinding=args.blind,
                                                          #rpf_params = params,
                                                          ),
                                                     )
                        
                        client.shutdown()
                        del cluster

                
                output['analysisCategories'] = label_map
                util.save(output, savefilename)
                print('saving', savefilename)

               



    elapsed = time.time() - tic
    printTime(elapsed)
    print(f"Events/s: {metrics['entries'] / elapsed:.0f}")
    



