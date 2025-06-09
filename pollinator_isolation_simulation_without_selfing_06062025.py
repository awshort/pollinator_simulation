"""
# Updated on June 09, 2025
# @authors: Carrie, Aidan, and John
# This code simulates pollinator isolation between Penstemon barbatus and virgatus without selfling
# Usage: python script.py <sm_val> <L> <QTL> <N> <gens> <reps> <sampled_gen> <seeds> <rr>
"""
import argparse
import random
import sys
import os

#################################################################

sm_val = float(sys.argv[1])
a1 = 1-sm_val
a2 = 1-a1
L = int(sys.argv[2]) # number of loci
QTL = int(sys.argv[3]) # number of floral trait loci
N = int(sys.argv[4]) # population size
gens = int(sys.argv[5]) # number of generations
reps = int(sys.argv[6]) # number of replicates
sampled_gen = int(sys.argv[7]) # sampling every XX generations
seeds = int(sys.argv[8]) # number of offspring produced per crossing event
rr = float(sys.argv[9]) # recombination rate between QTL and adjacent neutral loci

##################################################################

#function to create initial population with half bee adapted homozygotes and half hummingbird adapted homozygotes
def setup_pop(N, L):
    bee_adapted_geno = [[0 for i in range(L)],[0 for i in range(L)]] #Number of individuals fixed for bee alleles
    hb_adapted_geno = [[1 for i in range(L)],[1 for i in range(L)]] #Number of individuals fixed for hummingbird alleles
    pop = [] # generates a list of matrixes
    for i in range(int(N/2)):
        pop.append(bee_adapted_geno) # 
        pop.append(hb_adapted_geno) # 
    return(pop)

#function to select QTL positions
def generate_QTL_positions_2(L, QTL):
    my_list = list(range(L)) #All positions
    value = int(L/QTL) #number of neutral loci adjacent to QTL
    every_other_element = my_list[::value] #Select QTL positions with "value" number of neutral loci between them
    return(every_other_element)

#function to genotype individual at QTL 
def genotype_QTL(seed_parent,QTL_positions):
    my_list = [sum(col) for col in zip(*seed_parent)] 
    positions_to_sum = QTL_positions
    sum_of_elements = 0
    for position in positions_to_sum:
        if position < len(my_list):
            sum_of_elements += my_list[position]

    return(sum_of_elements)  

#function to calculate likelihood of bee visitation

def F_bee(x):
    global a1
    global a2
    return a1*((a1/a2)**(-x))

#function to calculate likelihood of hummingbird visitation

def F_hb(x):
    global a1
    global a2
    return a1*((a1/a2)**(-(1-x)))

#function to select which pollinator individual is visited by

def pollination(x):
    rand = random.random()
    if rand <= F_bee(x):
        pollinator = 'bee' # if less than rand visited by HB
    elif rand >= (1-F_hb(x)):
        pollinator = 'HB' # if less than rand visited by bee
    else:
        pollinator = 'no_visit' #'no_visit'
    return(pollinator)

#create gametes
def create_gamete(parent):
	gamete1 = parent[0]
	gamete2 = parent[1]
	gamete=[]
	for i in QTL_positions:
		if random.random() < 0.5:
			gameteA =  gamete1
			gameteB =  gamete2
		else:
			gameteA =  gamete2	
			gameteB =  gamete2	
		gamete.append(gameteA[i])
		for j in range(1,int(L/QTL)):
			if random.random() < rr:
				# Swap the alleles at this position
				gamete.append(gameteB[i+j])
			else:
				gamete.append(gameteA[i+j])          
	return(gamete)

#create zygote
import random
def create_zygote(mom, dad):
    egg = create_gamete(mom)
    sperm = create_gamete(dad)
    kid = [egg,sperm]
    return(kid)

# main
outfile = open('Z_tally.pollinator_simulation.sm_val_{4}.L{0}.QTL{6}.N{1}.gens{2}.reps{3}.rr{5}.offsrping_per_mating_pair{7}.without_selfing.txt'.format(L, N, gens, reps, sm_val, rr, QTL, seeds), 'w')
outfile.write('rep\tgen\tZ\ttally\tfraction\n')
outfile2 = open('QTL_tally.pollinator_simulation.sm_val_{4}.L{0}.QTL{6}.N{1}.gens{2}.reps{3}.rr{5}.offsrping_per_mating_pair{7}.without_selfing.txt'.format(L, N, gens, reps, sm_val, rr, QTL, seeds), 'w')
outfile2.write('rep\tgen\tQ\ttally\tfraction\n')

QTL_positions = generate_QTL_positions_2(L, QTL)
for r in range(reps):

    # starting population
    Z_tallies = [0 for i in range((2 * L) + 1)] # 11 possible categories (the number of individuals with 2*L hummingbird alleles starts with number of individuals with 0 hummingbird alleles ends with number of individuals with all hummingbird alleles)
    Z_tallies[0] += N / 2 # half pop homozygous for bee
    Z_tallies[-1] += N / 2 # half pop homozygous for hummingbird
    
    # starting population
    Q_tallies = [0 for i in range((2 * QTL) + 1)] # 11 possible categories (the number of individuals with 2*L hummingbird alleles starts with number of individuals with 0 hummingbird alleles ends with number of individuals with all hummingbird alleles)
    Q_tallies[0] += N / 2 # half pop homozygous for bee
    Q_tallies[-1] += N / 2 # half pop homozygous for hummingbird

    # set up the parent generation
    parent_pop = setup_pop(N, L)

    for j in range(gens):

        if j%sampled_gen == 0:
            tot_tallies = sum(Z_tallies)
            for z in range(len(Z_tallies)):
                outfile.write(str(r) + '\t' + str(j) + '\t' + str(z) + '\t' + str(Z_tallies[z]) + '\t' + str(float(Z_tallies[z])/(float(tot_tallies))) + '\n')
            tot_tallies = sum(Q_tallies)
            for q in range(len(Q_tallies)):
                outfile2.write(str(r) + '\t' + str(j) + '\t' + str(q) + '\t' + str(Q_tallies[q]) + '\t' + str(float(Q_tallies[q])/(float(tot_tallies))) + '\n')


        offspring_pop = []
        hb_visited = []
        bee_visited = []
        Z_tallies = [0 for i in range((2*L)+1)]
        Q_tallies = [0 for i in range((2*QTL)+1)]

        #for n in range(N):
        for n in range(len(parent_pop)):

            # for each individual in parent pop do this
            seed_parent = parent_pop[n]

            # visit by bee or hb (or skip) based on seed parent genotypic score
            Z = genotype_QTL(seed_parent,QTL_positions) # Z = no. HB alleles
            x = float(Z)/(2*QTL) # x = genotypic score
            pollen_parent = pollination(x)
            
            if pollen_parent == 'no_visit':
                
                pass

            elif pollen_parent == 'bee':
                    
                    bee_visited.append(seed_parent)
                    
            elif pollen_parent == 'HB':
                    
                    hb_visited.append(seed_parent)
                    
        for n in range(len(bee_visited)):
                                
            seed_parent = bee_visited[n]
            pollen_donor = random.choice(bee_visited)
            
            for i in range(seeds):
                
                # generate diploid genotypes based on seed and pollen parents and send to next generation
                offspring = create_zygote(seed_parent, pollen_donor)
                offspring_pop.append(offspring)
                
        for n in range(len(hb_visited)):
                                
            seed_parent = hb_visited[n]
            pollen_donor = random.choice(hb_visited)
            
            for i in range(seeds):
                
                # generate diploid genotypes based on seed and pollen parents and send to next generation
                offspring = create_zygote(seed_parent, pollen_donor)
                offspring_pop.append(offspring)
            
        #Remove offspring until len of offspring is equal to N          
        for n in range(len(offspring_pop)- N): 
            
            #draw a random offspring
            culled = random.choice(offspring_pop)
            offspring_pop.remove(culled)     
        
        # find offspring genotype category and save to tallies
        for i in range(len(offspring_pop)):
            
            offspring = offspring_pop[i]
            off_Z = sum([sum(col) for col in zip(*offspring)])
            Z_tallies[off_Z] += 1
            off_Q = genotype_QTL(offspring,QTL_positions) # Z = no. HB alleles
            Q_tallies[off_Q] += 1
            
        #Replace parent with offspring pop
        parent_pop = offspring_pop
        
tot_tallies = sum(Z_tallies)
for z in range(len(Z_tallies)):
    outfile.write(str(r) + '\t' + str(j) + '\t' + str(z) + '\t' + str(Z_tallies[z]) + '\t' + str(float(Z_tallies[z])/(float(tot_tallies))) + '\n')
tot_tallies = sum(Q_tallies)
for q in range(len(Q_tallies)):
    outfile2.write(str(r) + '\t' + str(j) + '\t' + str(q) + '\t' + str(Q_tallies[q]) + '\t' + str(float(Q_tallies[q])/(float(tot_tallies))) + '\n')

