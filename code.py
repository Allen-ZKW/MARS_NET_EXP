#call related module
import json
import networkx as nx
from SetCoverPy import setcover
from random import shuffle
from copy import deepcopy
from  numpy import *

#define functions
def undirect(reax_list):
    buffer = set()
    for r in reax_list:
        buffer.add(r.split('_')[0])
    buffer = list(buffer)    
    return buffer

def batchdel(tar_list,del_list):#functions to delete batch of data
    for i in sorted(del_list,reverse=True):
        del tar_list[i]
    return tar_list

def importdata(dirpath):#function to import related data from file 
    datapath = dirpath + '/data'
    with open(datapath + '/KEGG_compounds.json','r') as f1:
        kegg_com = json.load(f1)#imformation about metabolite: CID+NAME+MW
    
    with open(datapath + '/reaction(plus_found_by_ec)_in_organisms.json','r') as f2:
        rea_org = json.load(f2)#imformation about different species: NAME+ABB+TAX+ORTHOLOGY+REACTION+METABOLITES

    with open(datapath + '/allreactions.json','r') as f:
        kegg_reax = json.load(f)#imformation about reaction

    with open(datapath + '/Net550_reaction_direction.json','r') as f4:
        direction = json.load(f4)#imformation about the directions of reactions

    nutriGrup = {}
    with open(datapath + '/nutriGrup.txt','r') as f6:
        for le in f6:
            l = le.strip().split('\t')
            nutriGrup[l[0]]=l[-1]
    #imformation about each species's nutritional type

    envir = [] #551 TID in table S1
    with open(datapath + '/envir.txt') as f7:
        for le in f7:
            l = le.strip()
            envir.append(l)
    
    vip = [] #vip metabolite
    with open(datapath + '/vip.txt') as f8:
        for le in f8:
            l = le.strip()
            vip.append(l)
    vip = list(set(vip))
    
    sl = []#collect the reactions which product Ribose by decomposing
    with open(datapath + '/toRiboseRC.txt','r') as f9:
        for le in f9:
            l = le.strip()
            sl.append(l)

    return kegg_com,kegg_reax,rea_org,direction,nutriGrup,envir,sl,vip

def mkmatrix(kegg_com,kegg_reax):#function to create matrix which record imformation from reactions
    CID = list(kegg_com.keys())# get CID of each reactions
    RID = list(kegg_reax.keys())# get RID of each metabolites
    S_R = zeros((len(kegg_reax.keys()),len(kegg_com)))# init matrix to record reactant
    S_P = zeros((len(kegg_reax.keys()),len(kegg_com)))# init matrix to record resultant
    for key in kegg_reax.keys():
        x = RID.index(key)
        for i in kegg_reax[key]['reants'].keys():
            y = CID.index(i)
            S_R[x,y] = 1
        for j in kegg_reax[key]['produxs'].keys():
            y = CID.index(j)
            S_P[x,y] = 1
    return RID,CID,S_R,S_P,kegg_reax

def mkseed(meta_ori,CID):#function to create an array which record the information about seed
    seed = []
    for i in meta_ori:
        seed.append(CID.index(i))
    x=zeros((len(CID),1))
    x[seed]=1# 0: not seed; 1: seed
    return x

def CoRibose(kegg_reax,S_R,S_P,sl,kegg_com):#collection all reactions which product Ribose
    sll = []
    for r in sl:
        if r in kegg_reax.keys():
            for c in kegg_reax[r]['produxs']:
                for n in kegg_com[c]['name']:
                    if 'ribose' in n.lower() or 'glucose'in n.lower() or 'fructose'in n.lower():
                        sll.append(r+'_F')
            for c in kegg_reax[r]['reants']:
                for n in kegg_com[c]['name']:
                    if 'ribose' in n.lower() or 'glucose'in n.lower() or 'fructose'in n.lower():
                        sll.append(r+'_R')        
    toribose = list(set(sll))
    return toribose

def filiter_background(envir,direction,rea_org,toribose,CarFixReax,RID,nutriGrup,kegg_reax):#filiter the background
    background = {}
    for o in envir:
        background[o] = {}
        background[o]['reax'] = []
        for r in rea_org[o]['reaction']:
            if r in RID:
                if r in direction.keys():
                    if direction[r] == 'D':
                        background[o]['reax'].append(r+'_F')
                    elif direction[r] == 'F':
                        background[o]['reax'].append(r+'_F')
                else:
                    background[o]['reax'].append(r+'_F')

        for r in rea_org[o]['reaction']:
            if r in RID:
                if r in direction.keys():
                    if direction[r] == 'D':
                        background[o]['reax'].append(r+'_R')
                    elif direction[r] == 'R':
                        background[o]['reax'].append(r+'_R')
                else:
                    background[o]['reax'].append(r+'_R')
        # record all reactions and their direction in all species

        del_list = []
        for r in background[o]['reax']:
            reax = r.split('_')[0]
            d = r.split('_')[1]
            if d == 'F':
                reants = list(kegg_reax[reax]['reants'].keys())
            if d == 'R':
                reants = list(kegg_reax[reax]['produxs'].keys())
            if 'C00011' in reants or 'C00288' in reants or 'C00237' in reants:# only take reaction which fix carbon
                if 'auto' not in nutriGrup[o].lower():
                    del_list.append(background[o]['reax'].index(r))# carbon fix reaction cannot react in non autotrophic species
                elif r.split('_')[0] not in  CarFixReax:
                    del_list.append(background[o]['reax'].index(r))
            if r == 'R09503_R':# this reaction cannot happen in nature
                del_list.append(background[o]['reax'].index(r))
        background[o]['reax'] = batchdel(background[o]['reax'],del_list)
        
        background[o]['meta'] = set()
        del_list = []
        for r in background[o]['reax']:
            reax = r.split('_')[0]
            background[o]['meta'] = background[o]['meta'] | set(kegg_reax[reax]['reants'].keys()) | set(kegg_reax[reax]['produxs'].keys())
            if r in toribose:
                del_list.append(background[o]['reax'].index(r))
        background[o]['reax'] = batchdel(background[o]['reax'],del_list)
        background[o]['meta'] = list(background[o]['meta'])
    return background

def netEXP(R,P,x,b,reax,CID):#function to expand seed to network
# reference:Joshua Goldford
    k=sum(x,axis=0)[0]
    k0=0
    while k>k0:
        k0 = k
        s_d = transpose(dot(R,x))[0]
        y = array(s_d==b).astype(int)
        x_new = array(dot(y,P)>0).astype(int)
        for i in range(len(x_new)):
            if x_new[i]==1:
                x[i] = 1
        k = sum(x,axis = 0)[0]
    buffer = [k for k,t in enumerate(x) if t==1]
    meta = [j for i,j in enumerate(CID) if i in buffer]
    meta_n = k#count the number of metabolites after this loop
    buffer = [k for k,t in enumerate(y) if t==1]
    reax_ = [j for i,j in enumerate(reax) if i in buffer]
    reax_n = float(sum(y))#reaction exists in last step must exist in the next step 
    return meta,meta_n,reax_,reax_n

def netEXP_S(R,P,x,b,all_reax,CID):#function to expand seed to network step by step
#very similar to the obove function
    k=sum(x,axis=0)[0]
    k0=0
    meta = {}
    reax = {}
    num = 0
    buffer = [k for k,t in enumerate(x) if t==1]
    meta[num] = [j for i,j in enumerate(CID) if i in buffer]
    reax[num] = []
    while k>k0:
        num = num+1
        k0 = k
        s_d = transpose(dot(R,x))[0]
        y = array(s_d==b).astype(int)
        x_new = array(dot(y,P)>0).astype(int)
        for i in range(len(x_new)):
            if x_new[i]==1:
                x[i] = 1
        buffer = [k for k,t in enumerate(x) if t==1]
        meta[num] = [j for i,j in enumerate(CID) if i in buffer]
        buffer = [k for k,t in enumerate(y) if t==1]
        reax[num] = [j for i,j in enumerate(all_reax) if i in buffer]
        k = sum(x,axis = 0)[0]
    meta_n = k
    reax_n = float(sum(y))      
    return meta,meta_n,reax,reax_n

def single_organism(envir,rea_org,direction,RID,CID,S_R,S_P,toribose,nutriGrup,meta_ori,CarFixReax):# network for each species
    net_org = {}#init expand result
    for o in envir:
        x = mkseed(meta_ori,CID)
        reax = []
        ZX = []
        FX = []
        for r in rea_org[o]['reaction']:
            if r in RID:
                loc = RID.index(r)
                if r in direction.keys():
                    if direction[r] == 'D':
                        ZX.append(loc)
                        reax.append(r+'_F')
                    if direction[r] == 'F':
                        ZX.append(loc)
                        reax.append(r+'_F')
                else:
                    ZX.append(loc)
                    reax.append(r+'_F')
        for r in rea_org[o]['reaction']:
            if r in RID:
                loc = RID.index(r)
                if r in direction.keys():
                    if direction[r] == 'D':
                        FX.append(loc)
                        reax.append(r+'_R')
                    if direction[r] == 'R':
                        FX.append(loc)
                        reax.append(r+'_R')
                else:
                    FX.append(loc)
                    reax.append(r+'_R')
        R = vstack((S_R[ZX,:],S_P[FX,:]))
        P = vstack((S_P[ZX,:],S_R[FX,:]))

        del_list = []# init list to record delete locations
        for i in reax:
            if i in toribose:
                del_list.append(reax.index(i))
            elif R[reax.index(i),CID.index('C00011')]==1 or R[reax.index(i),CID.index('C00288')]==1 or R[reax.index(i),CID.index('C00237')]==1:
                if 'auto' not in nutriGrup[o].lower():
                    del_list.append(reax.index(i))
                elif i .split('_')[0] not in  CarFixReax:
                    del_list.append(reax.index(i))
                elif i == 'R09503_R':
                    del_list.append(reax.index(i))

        R = delete(R,del_list,axis=0)
        P = delete(P,del_list,axis=0)
        reax = batchdel(reax,del_list)
        b = sum(R,axis=1)
        net_org[o] = {'name':o,'meta':{},'meta_n':{},'reax':{},'reax_n':{}}
        net_org[o]['meta'],net_org[o]['meta_n'],net_org[o]['reax'],net_org[o]['reax_n'] = netEXP (R,P,x,b,reax,CID)
    return net_org

def multy_organism(envir,rea_org,direction,RID,CID,S_R,S_P,toribose,nutriGrup,meta_ori,CarFixReax):#all species network expand
    ZX = []
    FX = []
    all_reax = []
    for o in envir:
        for r in rea_org[o]['reaction']:
            if r in RID:
                loc = RID.index(r)
                if r in direction.keys():
                    if (direction[r] == 'D' or direction[r] == 'F') and r+'_F' not in all_reax:
                        ZX.append(loc)
                        all_reax.append(r+'_F')
                elif r+'_F' not in all_reax:
                    ZX.append(loc)
                    all_reax.append(r+'_F')
    for o in envir:
        for r in rea_org[o]['reaction']:
            if r in RID:
                loc = RID.index(r)
                if r in direction.keys():
                    if (direction[r] == 'D' or direction[r] == 'R') and r+'_R' not in all_reax:
                        FX.append(loc)
                        all_reax.append(r+'_R')
                elif r+'_R' not in all_reax:
                    FX.append(loc)
                    all_reax.append(r+'_R')

    R = vstack((S_R[ZX,:],S_P[FX,:]))
    P = vstack((S_P[ZX,:],S_R[FX,:]))
    del_list = []
    for i in all_reax:
        if i in toribose:
            del_list.append(all_reax.index(i))
        elif R[all_reax.index(i),CID.index('C00011')]==1 or R[all_reax.index(i),CID.index('C00288')]==1 or R[all_reax.index(i),CID.index('C00237')]==1:
            if i .split('_')[0] not in  CarFixReax:
                del_list.append(all_reax.index(i))
        elif i == 'R09503_R':
            del_list.append(all_reax.index(i))

    R = delete(R,del_list,axis=0)
    P = delete(P,del_list,axis=0)
    all_reax  = batchdel(all_reax,del_list)
    b = sum(R,axis=1)
    x = mkseed(meta_ori,CID)
    allt = {}
    allt = {'name':'all','meta':{},'meta_n':{},'reax':{},'reax_n':{}}
    allt['meta'],allt['meta_n'],allt['reax'],allt['reax_n'] = netEXP(R,P,x,b,all_reax,CID)
    x = mkseed(meta_ori,CID)
    allT = {}
    allT = {'name':'all','meta':{},'meta_n':{},'reax':{},'reax_n':{}}
    allT['meta'],allT['meta_n'],allT['reax'],allT['reax_n'] = netEXP_S(R,P,x,b,all_reax,CID)
    return allt,allT

def create_com(net_org,kegg_reax):#build relationship between metabolites and reactions
    reax_com = []
    for j in net_org['reax']:
        buffer = j.split('_')
        if buffer[1] == 'F':
            reagent = list(kegg_reax[buffer[0]]['reants'].keys())
            resultant = list(kegg_reax[buffer[0]]['produxs'].keys())
        else:
            reagent = list(kegg_reax[buffer[0]]['produxs'].keys())
            resultant = list(kegg_reax[buffer[0]]['reants'].keys())
        for b in resultant:
            for a in reagent:
                reax_com.append([a,b,j])
    return reax_com

def create_net(allT,reax_com,rea_org,envir,meta_ori,bioSeed):# use networkx to build network of reactions
    meta = allT['meta']
    reax = allT['reax']
    netGraph = []
    for l in reax_com:
        m1 = min([n for n in meta.keys() if l [0] in meta[n]])
        m2 = min([n for n in meta.keys() if l [1] in meta[n]])
        rn = min([n for n in reax.keys() if l [2] in reax[n]])
        ro = [] 
        for o in rea_org.keys():
            if o in envir and o.startswith('T'):
                if l[2].split('_')[0] in rea_org[o]['reaction']:
                    ro.append(o)
        netGraph.append([l[0],m1,l[1],m2,l[2],rn,len(ro),ro])

    netGraph_flt = []
    include_ori = []
    ori = []
    buffer = {}
    for l in netGraph:
        if l[1] == l[3]-1 and l[3] == l[5]:#find out reactions which really enlarge the network
            netGraph_flt.append(l)
            if l[0] in meta_ori:
                include_ori.append(l[0])
    
    ori = include_ori
    for l in netGraph:    
        if l[0] in meta_ori and l[0] not in include_ori:
            buffer[l[0]] = l
            include_ori.append(l[0])
        elif l[0] in meta_ori and l[0] not in ori and l[5] < buffer[l[0]][5]:
            buffer[l[0]] = l
    
    for i in buffer.keys():
        netGraph_flt.append(buffer[i])

    directed_G = nx.DiGraph()#init a directed graph
    for l in netGraph_flt:   
        directed_G.add_edge(l[0],l[2])#add nodes and edges to the graph
    return netGraph_flt,directed_G  

def vip_Graph(directed_G,vip,meta_ori,kegg_reax,netGraph_flt):# create a sub-graph which has all conditions for vip-produce
    target = set(vip)
    source = meta_ori
    subGraph = nx.DiGraph()
    subGraph.add_nodes_from(source)
    path = {}
    while target != set():
        require = set()
        reactions = set()
        for i in source:
            for j in target:
                n = list(nx.all_simple_paths(directed_G, source=i, target=j))#find path from seed towards vip
                if n != []:
                    path[i+'_'+j] = n
        for x in path.keys():
            for y in path[x]:
                for num in range(len(y)-1):
                    subGraph.add_edge(y[num],y[num+1])
        for p in subGraph.edges():
            reactions = reactions | set([q[4] for q in netGraph_flt if q[0] == p[0] and q[2] == p[1]])
        for i in reactions:
            buffer = i.split('_')
            if buffer[1] == 'F':
                require = require | set(kegg_reax[buffer[0]]['reants'].keys())
            if buffer[1] == 'R':
                require = require | set(kegg_reax[buffer[0]]['produxs'].keys())
        target = require - set(subGraph.nodes())
    return path

def path_com(path,netGraph_flt,meta_ori,kegg_reax):# requirment and produce of each paths
    path_index = {}
    for key in path.keys():
        if len(path[key]) == 1:
            path_index[key] = path[key][0]
        else:
            for i in range(len(path[key])):
                path_index[key+'_'+str(i+1)] = path[key][i]
    require = {}
    support = {}
    for index in path_index.keys():
        y = path_index[index]
        reactions = set()
        require[index] = set()
        support[index] = set()
        for num in range(len(y)-1):
            p = [y[num],y[num+1]]
            reactions = reactions | set([q[4] for q in netGraph_flt if q[0] == p[0] and q[2] == p[1]])
        for i in reactions:
            buffer = i.split('_')
            if buffer[1] == 'F':
                require[index] = require[index] | set(kegg_reax[buffer[0]]['reants'].keys())
                support[index] = support[index] | set(kegg_reax[buffer[0]]['produxs'].keys())
            if buffer[1] == 'R':
                require[index] = require[index] | set(kegg_reax[buffer[0]]['produxs'].keys())
                support[index] = support[index] | set(kegg_reax[buffer[0]]['reants'].keys())
    cross = {}
    for i in support:
        for j in support[i]:
            if j not in cross:
                cross[j] = set()
                cross[j].add(i)
            else:
                cross[j].add(i)
    return path_index,support,require,cross

def greedy_Graph(support,require,revip,path_index,cross,directed_G,meta_ori,netGraph_flt):
    greedy_G = nx.DiGraph()# init a directed graph
    greedy_G.add_nodes_from(meta_ori)# add seeds to the graph
    need = revip - set(greedy_G.nodes)
    reants = set()
    produxs= set(meta_ori)
    
    while need != set():
        score = {}
        for i in need:
            for j in cross[i]:
                if j not in score.keys():
                    score[j] = len(support[j] & need) - len(require[j] - produxs)#choose reaction which can produce more and require less        
        new = sorted(score.items(),key=lambda x:x[1],reverse=True)[0][0]
        reants = reants | require[new]
        produxs = produxs | support[new]
        y = path_index[new]
        for num in range(len(y)-1):
            greedy_G.add_edge(y[num],y[num+1])
        need = (revip | reants) - produxs
        for key in cross.keys():
            cross[key] = cross[key] - set([new])

    reactions = set()
    for p in greedy_G.edges:
        reactions = reactions | set([q[4] for q in netGraph_flt if q[0] == p[0] and q[2] == p[1]])
    expReaction = list(reactions)
    return expReaction

def setcover_problem(envir,background,expReaction):
    data = []
    for o in envir:
        io = []
        for n in expReaction:
            if n in background[o]['reax']:
                io.append(1)
            else:
                io.append(0)
        data.append(io)
    data = array(data, dtype=bool)#create cover matrix
    data = data.T
    cost = ones(len(envir), dtype=int)# all species' cost equal to 1
    g = setcover.SetCover(data, cost, maxiters= 20,subg_maxiters=200)
    solution, time_used = g.SolveSCP()
    result = list(g.s)
    leastOrganism = []
    n = -1
    for i in result:
        n +=1
        if i == True:
            leastOrganism.append(envir[result.index(i,n)])
    return leastOrganism

def part_organism(leastOrganism,rea_org,direction,RID,CID,S_R,S_P,toribose,nutriGrup,meta_ori,CarFixReax):
    # net-expand for given specied
    ZX = []
    FX = []
    part_reax = []
    for o in leastOrganism:
        for r in rea_org[o]['reaction']:
            if r in RID:
                loc = RID.index(r)
                if r in direction.keys():
                    if (direction[r] == 'D' or direction[r] == 'F') and r+'_F' not in part_reax:
                        ZX.append(loc)
                        part_reax.append(r+'_F')
                elif r+'_F' not in part_reax:
                    ZX.append(loc)
                    part_reax.append(r+'_F')
    for o in leastOrganism:
        for r in rea_org[o]['reaction']:
            if r in RID:
                loc = RID.index(r)
                if r in direction.keys():
                    if (direction[r] == 'D' or direction[r] == 'R') and r+'_R' not in part_reax:
                        FX.append(loc)
                        part_reax.append(r+'_R')
                elif r+'_R' not in part_reax:
                    FX.append(loc)
                    part_reax.append(r+'_R')
    R = vstack((S_R[ZX,:],S_P[FX,:]))
    P = vstack((S_P[ZX,:],S_R[FX,:]))
    del_list = []
    for i in part_reax:
        if i in toribose:
            del_list.append(part_reax.index(i))
        elif R[part_reax.index(i),CID.index('C00011')]==1 or R[part_reax.index(i),CID.index('C00288')]==1 or R[part_reax.index(i),CID.index('C00237')]==1:
            if i .split('_')[0] not in  CarFixReax:
                del_list.append(part_reax.index(i))
        elif i == 'R09503_R':
            del_list.append(part_reax.index(i))

    R = delete(R,del_list,axis=0)
    P = delete(P,del_list,axis=0)
    part_reax  = batchdel(part_reax,del_list)
    b = sum(R,axis=1)
    x = mkseed(meta_ori,CID)
    part = {}
    part = {'name':o,'meta':{},'meta_n':{},'reax':{},'reax_n':{}}
    part['meta'],part['meta_n'],part['reax'],part['reax_n'] = netEXP(R,P,x,b,part_reax,CID)
    return part

def caculate_rate(result,background):#caculate some imformation of the result
    for key in result.keys():
        result[key]['reax_rate'] = []
        result[key]['meta_rate'] = []
        result[key]['reax_num'] = []
        result[key]['meta_num'] = []
        for o in result[key]['org']:
            reax_num = len(set(undirect(background[o]['reax'])) & set(undirect(result[key]['reax'])))
            meta_num = len(set(background[o]['meta']) &  set(result[key]['meta']))
            reax_rate = reax_num / len(set(undirect(background[o]['reax'])))
            meta_rate = meta_num / len(set(background[o]['meta']))
            result[key]['reax_num'].append(reax_num)
            result[key]['meta_num'].append(meta_num)
            result[key]['reax_rate'].append(reax_rate)
            result[key]['meta_rate'].append(meta_rate)
    return result

def single_vip(net_org_s,vip,background):#caculate vip rate by single-result
    vip_rate = {}
    for key in net_org_s.keys():
        vip_rate[key] = len(vip & set(net_org_s[key]['meta']))/len(vip & set(background[key]['meta']))
    return vip_rate

def output_csv(target,dirpath,name):
    import csv
    with open(dirpath + '/result/' + name + ".csv",'w') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(["leastOrganism","reax_rate", "meta_rate","meta_num","reax_num","meta_n","reax_n"])
        for key in target.keys():
            writer.writerow([target[key]["org"],target[key]["reax_rate"],target[key]["meta_rate"],target[key]["meta_num"],target[key]['reax_num'],target[key]["meta_n"],target[key]["reax_n"]])

def main():
    meta_ori = ['C00001','C00009','C00011','C00014','C00059','C00237','C00244','C00283','C00533','C00697','C01353','C01438','C01455','C01485','C14818','C00028','C00030','C00002','C00003','C00004','C00010','C00024','C00205','C02061','C22151','C22150','C00667','C00662','C00139','C00138']
    bioSeed = ['C00028','C00030','C00002','C00003','C00004','C00010','C00024','C02061','C22151','C22150','C00667','C00662','C00139','C00138']
    CarFixReax = ['R00134','R00267','R00345','R01196','R01197','R07157','R10866','R00344','R00742','R01859','R00024','R00214','R00216','R00341','R00345','R07157','R10243']
    dirpath = 'D:/火星农业与健康/MAR_LEAST'
    kegg_com,kegg_reax,rea_org,direction,nutriGrup,envir,sl,vip = importdata(dirpath)
    key_meta = list(set(vip)|set(bioSeed))
    RID,CID,S_R,S_P,kegg_reax = mkmatrix(kegg_com,kegg_reax)
    toribose = CoRibose(kegg_reax,S_R,S_P,sl,kegg_com)
    background = filiter_background(envir,direction,rea_org,toribose,CarFixReax,RID,nutriGrup,kegg_reax)
    net_org_s = single_organism(envir,rea_org,direction,RID,CID,S_R,S_P,toribose,nutriGrup,meta_ori,CarFixReax)
    allt,allT = multy_organism(envir,rea_org,direction,RID,CID,S_R,S_P,toribose,nutriGrup,meta_ori,CarFixReax)
    reax_com = create_com(allt,kegg_reax)
    netGraph_flt,directed_G = create_net(allT,reax_com,rea_org,envir,meta_ori,bioSeed)
    path = vip_Graph(directed_G,key_meta,meta_ori,kegg_reax,netGraph_flt)
    rekey_meta = set(key_meta) & set(directed_G.nodes)
    path_index,support,require,cross = path_com(path,netGraph_flt,meta_ori,kegg_reax)
    expReaction = greedy_Graph(support,require,rekey_meta,path_index,cross,directed_G,meta_ori,netGraph_flt)
    vip_rate = single_vip(net_org_s,set(vip),rea_org)
    leastOrganism = []
    for i in range(1000):
        shuffle(envir)
        leastOrganism.append(setcover_problem(envir,background,expReaction))
    result = {}
    num = 0
    for i in leastOrganism:
        num += 1
        result[num] = part_organism(i,rea_org,direction,RID,CID,S_R,S_P,toribose,nutriGrup,meta_ori,CarFixReax)
        result[num]['org'] = i
    result = caculate_rate(result,background)
    output_csv(result,dirpath,'final')