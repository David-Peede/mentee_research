import sys
import numpy as np


# eval_accuracy measures the performance of an HMM against (Type 1) true introgressed sites
# INPUT:
#       true_windows = numpy array of windows (tiw) containing introgressed sites measured by percentage coverage.
#                      For purposes of clarity will ignore incomplete window edges as they result in Ns
#                      Treated within this method to only account for 100% introgressed windows
#       tested_HMM = numpy gamma matrix returned from HMM creation. Represents likelihood of introgression at position
#                    Should be pre-exponentiated (float probabilities not logs). Input is full numpy array
#       normalized = Boolean true/false determining behavior of algorithm. Default True
#                    False will HMM "guesses" as presented: are they over the threshold or not?
#                    True will adjust HMM "guesses", treating the highest-scoring window as 1 and normalizing other
#                    probabilities accordingly, then evaluating them in relation to the threshold.
#       threshold = float representing level of certainty HMM has of true state being introgressed.
#                   at or above threshold in the gamma matrix is evaluated as an "introgression guess" of the HMM
# OUTPUT: performance = an numpy array of floats [false_pos_r, miss_rate, sensitivity, specificity]
#       false_pos_r: "false positive" rate or (# false positives / (# false positives + # true negatives))
#                    probability true negative will test positive given the model (false alarm)
#       miss_rate:   "false negative" rate or (# false negatives / (# false negatives + # true positives))
#                    probability true positive will test negative given the model (expected often, use example)
#       sensitivity: "true positive" rate or (# true positives / (# true positives + # false negatives))
#                    probability true positive will test positive given the model
#       specificity: "true negative" rate or (# true negatives / (# true negatives + # false positives))
#                    probability true negative will test negative given the model (window edges treated as negative)
def eval_accuracy(true_windows, tested_HMM, normalized, threshold):
    tiw = true_windows  # matrix of windows containing introgressed sites (1. = 100% coverage)
    gamma = tested_HMM[:, 1]  # takes the second column of gamma, which has the Cs!
    # normalized = False  # set to True by default, comment this out for flexible use
    # threshold = .9  # set to .9 (90%) by default, comment out for flexible use
    # values preset to -1 to easily show errors
    false_pos_r = -1
    miss_rate = -1
    sensitivity = -1
    specificity = -1


    # IF we need to normalize, change values of gamma such that the highest-scoring site is set to 100%
    if normalized:
        # 1. Find the value of the highest-scoring site in gamma
        max_prob = np.amax(gamma)
        # 2. Multiply every probability in gamma (including the highest-scoring site) by 1/max_prob
        gamma = gamma * (1 / max_prob)

    # Compare gamma and tiw to see if they have the same number of elements (windows)
    if len(tiw) == len(gamma):

        fp = 0  # number of false positives
        fn = 0  # number of false negatives
        tp = 0  # number of true postivies
        tn = 0  # number of true negatives

        # loop through gamma and record the number of true/false positives/negatives
        for w in range(len(tiw)):
            # Underlying window is 100% introgressed
            if tiw[w] == 1.:
                # true positive
                if gamma[w] >= threshold:
                    tp += 1
                # false negative
                else:
                    fn += 1
            # Underlying window is not 100% introgressed
            elif 0 <= tiw[w] < 1.:
                # false positive
                if gamma[w] >= threshold:
                    fp += 1
                # true negative
                else:
                    tn += 1
            # something went wrong and tiw shows a value below zero or above 1
            else:
                print("ERROR in eval: window shows introgression percentage below zero or above 1")

        # check at the end to see if tp + tn + fp + fn == number of windows
        if tp + tn + fp + fn != len(tiw):
            print("ERROR in eval: true positives and true negatives do not sum to the number of windows")
            print(tp + tn + fp + fn)

        # Calculate rates for performance
        print("#tp: " + str(tp) + "\t#tn: " + str(tn), "\t#fp: " + str(fp), "\t#fn: " + str(fn))
        false_pos_r = fp / (fp + tn)
        miss_rate = fn / (fn + tp)
        sensitivity = tp / (tp + fn)
        specificity = tn / (tn + fp)

    else:
        print("ERROR IN EVAL: C column of gamma and true introgressed window array have different number of windows")

    performance = np.array([false_pos_r, miss_rate, sensitivity, specificity])
    return performance
