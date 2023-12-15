                                                                  Team 7 - Anomaly sound detection for predictive maintenance

We have used the dataset of Electrical engine availble at: https://www.idmt.fraunhofer.de/en/publications/datasets/isa-electric-engine.html and Made inference using Arduino Nano 33 BLE.

This Project has Three stages:

 Stage_1: As The device we are using is limited to sampling rate of 16000, we have downsampled the dataset for which we have used Resmple.m and as the dataset is small to increase the 
          efficiency of our design we have augmented it using reconstruction and noise addition (Which has been done using reconstruction.m).

 Stage_2: In this we have extracted the features from audio using a preprocessing unit with technique MFCC. For which the code is availble in c language format and arduino format in
          path: Project\MFCC. and input to it has given in .h format by converting .wav files into ,h files using header_creation.m file.

 Stage_3: In this the extracted features are given to classifer:
	  We have created two different classifiers for anomaly detection one with autoencoder design kind and the other with SVM based classification. Using Auto encoder we have found an 
	  accuracy of 69% and Using SVM was 93%. Codes for this are availble at Project\Python_Code

         (For this we have a Base Model created using Edge Impulse GUI, in which it is using MFCC as preprocessing unit and K-type Clustering method to classfy the input into anomaly or     
          normal. Code for this availble at Project\ei-haoml_project-arduino).

