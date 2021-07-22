import cobra, numpy, pandas, pickle, datetime
import multiprocessing, multiprocessing.pool

def manual_ko(working_gene):

    '''
    This function iterates genes, identifies reactions containing that particular gene, modifies boundaries of that reaction and then finally runs FBA to compute new objective value
    '''

    #! detect reations involved to that particular gene
    affected_reactions = []
    for reaction in model.reactions:
        local_reaction_genes = [gene.id for gene in reaction.genes]
        if working_gene in local_reaction_genes:
            affected_reactions.append(reaction)

    #! constrain reaction fluxes
    printt('before copy')
    #new_model = model.copy()
    #printt('after copy')
    #for working_reaction in affected_reactions:

    #    new_model.reactions.get_by_id(working_reaction.id).lower_bound = 0.
    #    new_model.reactions.get_by_id(working_reaction.id).upper_bound = 0.

    #! run FBA
    #new_solution = new_model.optimize()
    #print('objective value for {}: {}'.format(working_gene, new_solution.objective_value))

    #return new_solution.objective_value

def printt(message):

    print(datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S \t {}".format(message)))

    return None

###
### MAIN
###

###
### 1. user defined variables
###
model_file = '/Users/adrian/projects/riia/data/model/Recon3DModel_301.mat'
results_file = '/Users/adrian/projects/riia/results/manual_ko.pickle'
threads = 4

cobra_config = cobra.Configuration()
cobra_config.processes = threads
print(cobra.Configuration())

###
### 2. reading model
###
printt('reading model')
model = cobra.io.load_matlab_model(model_file)
printt('reading done')

### 3. running manual KOs
gene_ids = [gene.id for gene in model.genes]
tasks = gene_ids[:12]
print(tasks)

for task in tasks:
    manual_ko(task)

#hydra = multiprocessing.pool.Pool(threads)
#output = hydra.map(manual_ko, tasks)
#print(output)
