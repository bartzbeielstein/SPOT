options(spot.run.full.test = F)

if(!getOption("spot.run.full.test")){
    skip("Skipping shorter Tests due to option spot.run.full.tests = F in test.000TestSetup.R")
}