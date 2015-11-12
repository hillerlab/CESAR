import os  #for handling directories
import time
import params
import helpers
import model

def _find_smallest_codon_sub_prob():
    #function that finds the codonpair with the smallest sub prob
    #THIS IS NOT NEEDED ANY WHERE IN BUILDING THE MODEL. 

    
    count = 0
    smallest_value = 1
    codon_pair = ("-", "-")
    for k, v in params.codon_subs_probs.items():
        print(k, len(v))
        for k2, v in params.codon_subs_probs[k].items():
            if v < smallest_value:
                smallest_value = v
                codon_pair = (k,k2)
    return (smallest_value, codon_pair[0], codon_pair[1])


def _find_edge_codon_sub_prob(s = {"smallest", "largest"}):
    #function that finds the smallest/largest substitution prob
    #THIS IS NOT NEEDED ANY WHERE IN BUILDING THE MODEL. 
    count = 0
    if s == "smallest": edge_value = 1
    if s == "largest": edge_value = 0 
    codon_list = []
    for k in params.codon_subs_probs.keys():
        for k2, v in params.codon_subs_probs[k].items():
            count +=1
            if v == edge_value:
                codon_list.append((k,k2))
            if s=="smallest" and v < edge_value:
                edge_value = v
                codon_list = []
                codon_list.append((k,k2))
            if s=="largest" and v > edge_value:
                print("new largest value found" + str(v))
                edge_value = v
                codon_list = []
                codon_list.append((k,k2))
    print("number of pairs examined: " + str(count))
    return (edge_value, codon_list)


def _codon_sub_prob_avg_median(percentile):
    #THIS IS NOT NEEDED ANY WHERE IN BUILDING THE MODEL. 
    values_list = []
    for k in params.codon_subs_probs.keys():
        for k2, v in params.codon_subs_probs[k].items():
            values_list.append(v)
    values_list.sort()
    avg = sum(values_list)/len(values_list)
    median =values_list[len(values_list)/2]
    percentile = values_list[int(len(values_list)*percentile)]
    return (avg, median, percentile)


def aminoacid_edge_sub_prob(s = {"smallest", "largest"}):
    #THIS IS NOT NEEDED ANY WHERE IN BUILDING THE MODEL. 
    count = 0
    if s == "smallest": edge_value = 1
    if s == "largest": edge_value = 0 
    ac_list = []
    for k in params.aminoacid_subs_probs.keys():
        for k2, v in params.aminoacid_subs_probs[k].items():
            count +=1
            if v == edge_value:
                ac_list.append((k,k2))
            if s=="smallest" and v < edge_value:
                edge_value = v
                ac_list = [(k,k2)]
            if s=="largest" and v> edge_value:
                edge_value = v
                ac_list = [(k,k2)]
    print("number of pairs examined: " + str(count))
    return (edge_value, ac_list)

def _aminoacid_sub_prob_avg_median(percentile):
    #THIS IS NOT NEEDED ANY WHERE IN BUILDING THE MODEL. 
    values_list = []
    for k in params.aminoacid_subs_probs.keys():
        for k2, v in params.aminoacid_subs_probs[k].items():
            if k==k2: continue
            values_list.append(v)
    values_list.sort()
    print(values_list)
    print("length" + str(len(values_list)))
    avg = sum(values_list)/len(values_list)
    median =values_list[len(values_list)/2]
    percentile = values_list[int(len(values_list)*percentile)]
    return (avg, median, percentile)




def compare_two_sequences(s1, s2):

    #THIS IS NOT NEEDED ANY WHERE IN BUILDING THE MODEL. 
    if len(s1) != len(s2):
        print("the sequences are not of equal length")
    diff1 = ""
    diff2 = ""
    location = 0
    in_diff_zone = False
    for i in range(len(s1)):
         if s1[i] != s2[i]:
             if not in_diff_zone:  #entering difference zone now
                 in_diff_zone = True
                 location = i
                 diff1 = ""
                 diff2 = ""
             diff1 += s1[i]
             diff2 += s2[i]
         
         else:  #symbols match
             if in_diff_zone:
                 #if true then i is the first symbol after a exiting
                 #a difference area
                 in_diff_zone = False
                 print("Difference at location " + str(location))
                 print(diff1)
                 print(diff2)
    if in_diff_zone:
        print("Difference at location " + str(location))
        print(diff1)
        print(diff2)




def locate_frameshifts(queryseq, path):
    fs = []
    for i in range(len(path)-1):
        if  helpers.is_frameshift(path, i):
            location_in_seq = helpers.map_state_to_query_seq_position(path, i)
            fs.append((location_in_seq, path[i], path[i+1]))
        if helpers.emits_stop_codon(queryseq, path, i):
            location_in_seq = helpers.map_state_to_query_seq_position(path, i)
            fs.append((location_in_seq, path[i], "stopcodon"))
    return fs



def run_experiment(folders, idnum, output_name, paramstuple):
    '''
    This is the top level function that handles the experiment.

    It takes as arguments:
        folders: list of folders
        idnum: an identifier used to distinguish a experiment from another
            the idnum is used in naming the output file name, and
            in results table.
        output_name: name of output file where the alignments (and other
                results are written)
        paramstuple: the tuple of fs_prob, ci_prob, ci2_prob, total_cd_prob,
                substitution matrix
    '''
    
    #print("order of folders here is : " + str(folders))
    f = file(output_name, 'wr+')
    _write_params(f, paramstuple)
    stats = {}
    ratios = {}
    aggregated_scores = {}
    for folder in folders:
        #test_batch returns a tuple of 3 objects
        stats[folder], aggregated_scores[folder], ratios[folder] = test_batch(folder, output_name, paramstuple)
    _write_stats(output_name, stats)
    _write_aggregated_scores(output_name, aggregated_scores)
    _write_ratios(output_name, ratios)
    _type_averages_to_int(aggregated_scores)
    _write_results_oneline(output_name, idnum, paramstuple,
                           stats, aggregated_scores)
    f.close()
    return stats, aggregated_scores



def _write_stats(filename, stats):
    #function that writes statistics about a
    #folder into an experiment's output file
    f = file(filename, 'a')
    for folder in stats.keys():  
        f.write("\nSTATS for folder \"" + folder + "\":")
        f.write("\ncorrect alignment predictions: " +
                str(stats[folder][6]))
        f.write("\nwrong alignment predictions: " +
                str(stats[folder][7]))
        f.write("\nfiles w/ mixed results in recognizing the exon: " +
                str(stats[folder][8]))
                
        f.write("\n0 frameshifts: " + str(stats[folder][0]))
        f.write("\n1 frameshift: " + str(stats[folder][1]))
        f.write("\n2 frameshifts: " + str(stats[folder][2]))
        f.write("\n3 frameshifts: " + str(stats[folder][3]))
        f.write("\n4 or more frameshifts: " + str(stats[folder][4]))
        f.write("\ninstances skipped: "+ str(stats[folder][5]))
        f.write("\n..................................\n\n")
    f.flush()


def _write_aggregated_scores(filename, scores):
    f = file(filename, 'a')
    for folder in scores.keys():  
        f.write("\nALIGNMENT SCORES (prob of viterbi/prob of original alignment) "+
                "for folder \"" + folder + "\":")
        f.write("\ninstances with alignment score 1: " +
                str(scores[folder][0]))
        f.write("\ninstances with alignment score between 1 and 10: " +
                str(scores[folder][1]))
        f.write("\ninstances with alignment scores between 10 and 100: " +
                str(scores[folder][2]))
        f.write("\ninstances with alignment scores over 100: "+
                str(scores[folder][3]))
        f.write("\narithmatic average of the scores: "+
                str(scores[folder][4]))
        f.write("\nlog of the arithmatic average of the scores: "+
                str(scores[folder][5]))
        f.write("\nharmonic mean of the scores: "+
                str(scores[folder][6]))
        f.write("\n..................................\n\n")
    f.flush()


def _write_ratios(filename, ratios):
    #the ratios here refer to prob(viterbi-alignment)/prob(original-alignment)
    #this function writes the ratio of every instance in the experiment folder.
    f = file(filename, 'a')
    for folder in ratios.keys():
        ratios_float = [float(x) for x in ratios[folder]]
        ratios_float.sort(reverse= True)
        f.write("\nRATIOS FOR FOLDER \"" + folder + "\":\n")
        f.write(str(ratios_float))
    f.write("\n..................................\n\n")
    

def _type_averages_to_int(scores):
    #function to get rid of long floating point numbers
    for folder in scores.keys():
        scores[folder][5] = int(scores[folder][5])
    return scores

def _write_results_oneline(filename, idnum, paramstuple, stats, scores):
    #function that writs a summary of the results of the experiment as one line
    #to an excel-like file.
    f = file(filename, 'a')
    f.write(str(idnum) + ",")
    f.write(str(paramstuple[0:5] + (paramstuple[-1],)) + ",")
    for folder in stats.keys():
        f.write(str(stats[folder]) + ",")
        f.write(str(scores[folder]) + ",")
    f.write("\n")




def _write_query_details(query_index, f, instance, result,
                         frame_shifts, alignments_match,
                         score_original_alignment, score_ratio, 
                         path_decoded_from_original,
                         acc_site_correct, do_site_correct):
        #function that writes details about aligned pair into the output file
        queries_names = instance[6]
        f.write("***QUERY " + str(queries_names[query_index]))
        f.write("\noriginal alignment:")
        offset = (helpers.calculate_offset(instance[1][query_index]) -
                  helpers.leading_intron_in_sequence(instance[0]))
        f.write("\n" + (" "* offset) + instance[0])
        f.write("\n" + instance[1][query_index])
        f.write("\npredicted alignment:\n")
        f.write(result[1][0] + "\n")
        f.write(result[1][1] + "\n")
        f.write("alignments match: " +
                str(alignments_match))
        f.write("\nACC site recognized correctly: " + str(acc_site_correct))
        f.write("\nDO site recognized correctly:  " + str(do_site_correct))
        f.write("\nframeshifts count: " + str(len(frame_shifts)))
        f.write("\nframeshifts: " + str(frame_shifts))
        f.write("\nViterbi path probability: " + str(result[2]))
        f.write("\ncompile time: " + str(result[3]))
        f.write("\nviterbi calculation time: " +str(result[4]))

        #f.write("\nscore of original alignment: " + str(score_original_alignment))
        #f.write("\nscore ratio: "+ str(score_ratio))
        #f.write("\nscore ratio class: " + str(helpers.get_ratio_class(score_ratio)))
        #if result[0] == path_decoded_from_original:
         #   f.write("\nPAIR_TO_PATH FUNCTION WORKS HERE")
            #f.write("viterbi path: " + str(result[0]))
        #else:
         #   f.write("\nFAIL: PAIR_TO_PATH FUNCTION DOESN'T WORK HERE")
          #  f.write("\nviterbi path: " + str(result[0]))
           # f.write("\ndecoded path: " + str(path_decoded_from_original))
           # difference = []
            #for i,s in enumerate(result[0]):
            #    if i>= len(path_decoded_from_original) or path_decoded_from_original[i] !=s:
            #        difference.append((i, s))
            #f.write("\n difference: " + str(difference))
        f.write("\n\n")     



def test_batch(folder_path, output_name, paramstuple):
    if folder_path[-1] != "/":
        folder_path += "/"
    stats = [0] * 11
    score_ratios = []
    f = file(output_name, 'a')
    f.write("--------------------------------")
    f.write("\n\nCurrently processing folder " + folder_path + "\n\n")
    file_names = os.listdir(folder_path)
    for filename in file_names:
        instance = helpers.parse_instance(folder_path + filename)
        print("Currently testing: "+ filename)
        print("length of reference sequence: " +  str(len(instance[2])))
        print("number of query sequences: " +  str(len(instance[3])))
        #modelname = filename[0:len(filename)-4]
        modelname  = filename
        results, m = test_instance(modelname, instance, paramstuple)
        f.write("\n\n\n" + filename + ":\n")
        match_results = []
        for query_index, result in enumerate(results):
            frame_shifts = locate_frameshifts(instance[3][query_index],
                                              result[0])
            if frame_shifts == []:
                stats[0] += 1
            elif len(frame_shifts) == 1:
                stats[1] += 1
            elif len(frame_shifts) == 2:
                stats[2] += 1
            elif len(frame_shifts) == 3:
                stats[3] += 1
            else:  #4 or more frame-shifts
                stats[4] += 1
            '''
            print("\nVITERBI PATH: "+ str(result[0]))
            print("\nPREDICTED ALIGNMENT:")
            print("\n" + result[1][0])
            print("\n" + result[1][1])
            print("\nORIGINAL ALIGNMENT:")
            print("\n" + instance[0])
            print("\n" + instance[1][query_index])
            '''
            alignments_match = _alignments_match((result[1][0], result[1][1]),
                                             (instance[0], instance[1][query_index]))

            match_results.append(alignments_match)
            if alignments_match:
                stats[6] += 1
            else:
                stats[7] += 1
           
            #print("in analytics: result[1] " + str(result[1]))
            path_decoded_from_original = helpers.pair_to_path(m, result[1])
            score_original_alignment = 0
            #( model.score_alignment_against_model(m,(instance[0],
            #                                          instance[1][query_index]))
            if score_original_alignment == 0:
                score_ratio = 0
            else:
                score_ratio = result[2]/score_original_alignment


            score_ratios.append(score_ratio)

            acc_site_index_original = helpers.get_exon_start(instance[1][query_index])
            acc_site_index_predicted = helpers.get_exon_start(result[1][1])
            do_site_index_original = helpers.get_exon_end(instance[1][query_index])
            do_site_index_predicted = helpers.get_exon_end(result[1][1])
            acc_site_correct = (acc_site_index_original ==
                                acc_site_index_predicted)
            do_site_correct = (do_site_index_original ==
                               do_site_index_predicted)
            if acc_site_correct:
                stats[9] += 1
            if do_site_correct:
                stats[10] += 1
                
            _write_query_details(query_index, f, instance, result,
                                 frame_shifts, alignments_match,
                                 score_original_alignment, score_ratio,
                                 path_decoded_from_original,
                                 acc_site_correct,
                                 do_site_correct)
        f.write("\n\nmixed results: ")
        if sum(match_results) == 0 or sum(match_results) == len(match_results):
            f.write("False")
        else:
            f.write("True")
            stats[8] += 1

    f.flush()
    aggregated_ratios = helpers.aggregate_ratios(score_ratios)
    return stats, aggregated_ratios, score_ratios



def _alignments_match(pair1, pair2):
    #function is used to compare the original alignment to the predicted alignment
    if (pair1[0].strip() != pair2[0].strip() or
        pair1[1].strip() != pair2[1].strip()):
        return False
    else:
        return True
    
    

def generate_all_param_combinations():
    ''' function used to generate cluster commands file with
        all the combinations of parameters that have to be tested
    '''
    filename = "params_combinations.txt"
    fs_values = [0.001]
    ci_values = [0.001, 0.01, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3]
    ci2_values = [0.001, 0.01, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.5, 0.55]
    cd_values =  [0.001, 0.01, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3]
    nti_nti_values = [0.25]

    f = file(filename, "wr+")
    counter = 0
    for fs_p in fs_values:
        for ci_p in ci_values:
            for ci2_p in ci2_values:
                for cd_p in cd_values:
                    for nti_p in nti_nti_values:
                        f.write(str(counter) + "a")
                        f.write(":(")
                        f.write(str(fs_p) + "," + str(ci_p) + "," +
                                str(ci2_p) + "," + str(cd_p) + "," +
                                str(nti_p) + ", 'BLOSUM'")
                        f.write(")\n")
                        f.write(str(counter) + "b")
                        f.write(":(")
                        f.write(str(fs_p) + "," + str(ci_p) + "," +
                                str(ci2_p) + "," + str(cd_p) + "," +
                                str(nti_p) + ", 'ETH'")
                        f.write(")\n")
                        counter += 1
    f.close()



def generate_all_param_combinations2():
    
    ''' function used to generate cluster commands file with
        all the combinations of parameters that have to be tested
        Another version
    '''

    filename = "params_combinations.txt"
    fs_values = [0.0000001, 0.000001, 0.00001, 0.0001, 0.001, 0.01, 0.1]
    ci_values =  [0.1]
    ci2_values = [0.8]
    cd_values =  [0.01]
    nti_nti_values = [0.25]
    multiple_cd_factors = [0.432, 0.276, 0.208, 0.164, 0.147, 0.133, 0.127, 0.123, 0.118]

    f = file(filename, "wr+")
    counter = 0
    for fs_p in fs_values:
        for ci_p in ci_values:
            for ci2_p in ci2_values:
                for cd_p in cd_values:
                    for nti_p in nti_nti_values:
                        '''
                        f.write(str(counter) + "-1")
                        f.write(":(")
                        f.write(str(fs_p) + "," + str(ci_p) + "," +
                                str(ci2_p) + "," + str(nti_p) + "," +
                                str(cd_p) + ", 'BLOSUM'")
                        f.write(")\n")
                        
                        f.write(str(counter) + "-2")
                        f.write(":(")
                        f.write(str(fs_p) + "," + str(ci_p) + "," +
                                str(ci2_p) + "," + str(nti_p) + "," +
                                str(cd_p) + ", 'ETH'")
                        f.write(")\n")
                        '''


                        multiple_cd_probs = ""
                        for v in multiple_cd_factors:
                            multiple_cd_probs += str(v*cd_p) + ","
                        

                        f.write(str(counter) + "a")
                        f.write(":(")
                        f.write(str(fs_p) + "," + str(ci_p) + "," +
                                str(ci2_p) + "," + str(nti_p) + "," +
                                str(cd_p) + "," + multiple_cd_probs +
                                "'BLOSUM'")
                        f.write(")\n")

                        f.write(str(counter) + "b")
                        f.write(":(")
                        f.write(str(fs_p) + "," + str(ci_p) + "," +
                                str(ci2_p) + "," + str(nti_p) + "," +
                                str(cd_p) + "," + multiple_cd_probs +
                                "'ETH'")
                        f.write(")\n")

                        
                        counter += 1
    f.close()



def generate_cluster_file(params_file, folder1, folder2, id_prefix):
    
    ''' function used to generate cluster commands file with
        all the combinations of parameters that have to be tested
    '''

    f = file(params_file)
    output = file("3Nov_cluster_commands.txt", 'wr+')
    lines = f.readlines()
    for l in lines:
        parts = l.split(":")
        idnum = parts[0]
        p = parts[1].strip()
        output.write("python CESAR/run_experiment.py " +
                     folder1 + " " + folder2 + " " + id_prefix + idnum +
                     " \"" + p + "\"\n")
    output.close()



def _write_params(f, paramstuple):
    #helper function used to write the parameters of an experiment to the
    #output file of that experiment
    f.write("fs_prob = " + str(paramstuple[0]) + "\n")
    f.write("ci_prob = " + str(paramstuple[1]) + "\n")
    f.write("ci2_prob = " + str(paramstuple[2]) + "\n")
    f.write("total_cd_probs = " + str(paramstuple[3]) + "\n")
    f.write("codon substitution matrix: " + paramstuple[4] + "\n")
    f.write("\n\n")
    f.flush()
    




def test_instance(modelname, instance, paramsconfig):
    '''
    THis function builds the model and predicts the alignment for
    one instance. AN instance here means 1 ref sequence with 1 OR MORE query
    sequences.
    '''
    p = params.Params_config()
    p.set_params(paramsconfig)
    p.set_sequence(instance[2], instance[4], instance[5])
    m = model.HMM(modelname, p)
    results = []
    for query_seq in instance[3]:
        #print("current query sequence: " + query_seq)
        results.append(m.viterbi(query_seq, False))
    return results, m

    



def test_instance_from_file(filename, paramsconfig):
    '''
    THis function parses an instance from a file and then builds the model
    and predicts the alignment for
    one instance.

    AN instance here means 1 ref sequence with 1 or more query
    sequences
    '''
    instance = helpers.parse_instance(filename)
    modelname = filename[0:len(filename)-4]
    t1 = time.time()
    test_instance(modelname, instance, paramsconfig)
    t2 =time.time()
    print("total processing time: " + str(t2-t1) + " seconds.")
    
        
                
                
