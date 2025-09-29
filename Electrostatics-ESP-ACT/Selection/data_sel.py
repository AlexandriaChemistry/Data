#!/usr/bin/env python3
import random
# generate random test and train script
data = [
    "water#water|Train",
    "water#bromide|Train",
    "water#chloride|Train",
    "water#fluoride|Train",
    "water#lithium|Train",
    "water#sodium|Train",
    "water#potassium|Train",
    "acetate#acetate|Train",
    "acetate#lithium|Train",
    "acetate#sodium|Train",
    "acetate#potassium|Train",
    "acetate#water|Train",
    "acetate#fluoride|Train",
    "acetate#chloride|Train",
    "acetate#bromide|Train",
    "methylammonium#methylammonium|Train",
    "methylammonium#fluoride|Train",
    "methylammonium#chloride|Train",
    "methylammonium#bromide|Train",
    "methylammonium#formate|Train",
    "methylammonium#acetate|Train",
    "methylammonium#water|Train",
    "methylammonium#lithium|Train",
    "methylammonium#sodium|Train",
    "ammonium#acetate|Train",
    "ammonium#ammonium|Train",
    "ammonium#bromide|Train",
    "ammonium#chloride|Train",
    "ammonium#fluoride|Train",
    "ammonium#water|Train",
    "ammonium#formate|Train",
    "ammonium#lithium|Train",
    "ammonium#sodium|Train",
    "ethylammonium#acetate|Train",
    "ethylammonium#bromide|Train",
    "ethylammonium#chloride|Train",
    "ethylammonium#ethylammonium|Train",
    "ethylammonium#fluoride|Train",
    "ethylammonium#water|Train",
    "ethylammonium#formate|Train",
    "formate#formate|Train",
    "formate#lithium|Train",
    "formate#potassium|Train",
    "formate#sodium|Train",
    "formate#water|Train",
    "lithium#litium|Train",
    "sodium#sodium|Train",
    "potassium#potassium|Train",
    "fluoride#fluoride|Train",
    "chloride#chloride|Train",
    "bromide#bromide|Train",
    "sodium#chloride|Train",
    "lithium#fluoride|Train",
    "potassium#bromide|Train",
    "fluoride#chloride|Test",
    "sodium#potassium|Test",
    "ammonium#potassium|Test",
    "ethylammonium#lithium|Test",
    "ethylammonium#sodium|Test",
    "ethylammonium#potassium|Test",
    "formate#fluoride|Test",
    "formate#chloride|Test",
    "formate#bromide|Test",
    "methylammonium#potassium|Test",
    "potassium#chloride|Test",
    "potassium#fluoride|Test",
    "sodium#fluoride|Test",
    "sodium#bromide|Test",
    "lithium#chloride|Test",
    "lithium#bromide|Test",
    "lithium#potassium|Test",
    "lithium#sodium|Test",
    "fluoride#bromide|Test",
    "chloride#bromide|Test",
    "ammonium#ethylammonium|Test",
    "formate#acetate|Test",
    "ethylammonium#methylammonium|Test",
    "ammonium#methylammonium|Test"
]

cleaned_data = [entry.replace("|Train", "").replace("|Test", "") for entry in data]
train_data = [entry for entry in cleaned_data if "Train" in data[cleaned_data.index(entry)]]
test_data = [entry for entry in cleaned_data if "Test" in data[cleaned_data.index(entry)]]

combined_data = train_data + test_data
random.shuffle(combined_data)

test_count = int(0.30 * len(combined_data))

new_test_data = combined_data[:test_count]
new_train_data = combined_data[test_count:]
new_train_data = [entry + "|Train" for entry in new_train_data]
new_test_data = [entry + "|Test" for entry in new_test_data]

print("New Train Data:")
for entry in new_train_data:
    print(entry)

print("\nNew Test Data:")
for entry in new_test_data:
    print(entry)

