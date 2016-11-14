import re
import itertools


def read_file(file):
    rules = []
    with open(file) as rules_file:
        rules = rules_file.readlines()
    return rules


def split_by(text, delimiter):
    words = []
    for line in text:
        if not isinstance(line,list):
            words.append(line.split(delimiter))
        else:
            for x in line:
                words.append(x.split(delimiter))
    return words


def clean_rules(text):
    for row in text:
        for column in row:
            if column == '':
                row.remove(column)
            if column == '\n':
                row.remove(column)
    for row in text:
        for column in row:
            if column:
                if not column[0] == '{':
                    row.remove(column)
            else:
                row.remove(column)
    text = [row for row in text if row != []]
    return text


def make_rules(rules):
    words = []
    for row in rules:
        for line in row:
            line = line.replace('{', '')
            line = line.replace('}', '')
            line = line.replace(' ', '')
            line = line.replace('=>', ',')
            words.append(line)
    return words


def rules():
    words = read_file("newrules_MF_thresh_7_sort_lift.txt")
    words = split_by(words, '\t')
    words = split_by(words, '"')
    words = clean_rules(words)
    words = [x for x in words if x != []]
    words = make_rules(words)
    ruleset = split_by(words, ",")

    for rules in ruleset:
        rules[0], rules[len(rules)-1] = rules[len(rules)-1], rules[0]
    return ruleset
    # print(ruleset)


def clean_chosen(chosen):
    w = "\"x\"\n"
    spl = [list(y) for x, y in itertools.groupby(chosen, lambda z: z == w) if not x]
    chosen = []
    for x in spl:
        chosen.append([y[5:-2] for y in x])
    return chosen

def chosen():
    words = read_file("code/chosen.csv")
    words = clean_chosen(words)
    return words


def check(ruleset, chosen):
    correct = 0
    incorrect = 0
    for x in chosen:
        for y in ruleset:
            if set(y[1:]) <= set(x):
                if y[0] in x:
                    correct += 1
                else:
                    incorrect += 1
                    count = 0
                    for q in y:
                        if q in x:
                            count +=1
                    print(count/len(y))
    print('correct', correct)
    print('incorrect', incorrect)


def main():
    ruleset = rules()
    choices = chosen()
    check(ruleset, choices)
    #print(ruleset)
    #print('\n ------------------ \n')
    #print(choices)



main()