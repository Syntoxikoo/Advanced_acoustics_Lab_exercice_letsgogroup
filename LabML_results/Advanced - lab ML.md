## Conclusion drawn from Experiment
This labs presented two differents method in order to estimate the direction of arrival of a certain sound in a room within an angular resolution of 10°. The first method was a Convolutionnal Neural Network which is able to predict pretty accurately the DOA from a 10 channel STFT of the Intensity estimation. The second method (presented as a baseline to assess the performance of the CNN) is called SRP-PHAT which, basicaly, uses the signal received by the microphone two and the covariance between them to find the DOA.

The CNN shows great results while predicting the DOA from these STFTs of Intensity estimate, but several step was needed in order to be able to predict the DOA. Firstly simulating a lot of  room impulse response from  

Overall the CNN shows much more accurate prediction, and, once it was trained, was really fast to predict the DOA of a several signals. In comparison,  SRP HAT showed some large inacuracy while prediciting certain instance of the test set. In addition, has this model need to compute a lot of different quantity, it was quite slow for predicting a set of datapoint. As DOA is often use for real time application, this prediction can gain for more consistency in the results and more speed for this prediction. 
Nevertheless, this method has some downside, mainly that it needs to be trained on a lot of data, while the other algoritme can perform prediction from scratch. It is complicated to obtain quality data in quantity to train such model, in the present case we needed to simulate those data through room simulation using superposition principle, image source method and raytracing. Even if this simulation is to an extent close to the reality, it does not completly capturate the variance and uncertainty that we could find in real world measurement which could lead to big issue while asking the model to predict the DOA from real microphones datas

Training is computationnaly more expensive, weights take some memory usage even if in the presente case it is moderate, if the model would be wanted to perfom better and to have also finner resolution between angles, as well as being better in situation with a lot of reverberation
## Performance of the CNN model against SRP-PHAT estimation
The estimation of DOA using a neural network an intensity estimate provided really accurate results, with a high rate of accurate classification, when the network was misclassifying the DOA of a certain signal, the predicted one was really close to the true DOA, within 10° (1 resolution point) as we can observe in the figure above showing the true estimate and predicted one by the CNN. In the other hand, in the case of SRP-PHAT, while most of the prediction was close to the true DOA, for certain signal the results were completely off. It is hard to find the correct explanation for those misclassification, but as the signal received by the microphones is "white noise" convolved with the impulse response of the room, it can be due to some artifact while correlating this noise received by the differents microphones at differents time stamp.

In order to quantify by how much the CNN was better at predicting the DOA, the RMS error between true value and predicted one was computed for both model and plotted in the last figures. Results shows that the CNN has an error more than 10 times smaller ($0.057$ rad) than the SRP-PHAT algorithm ($0.727$ rad). This improvement is considerable regarding the fact that the second algorithm is a widely used algorithm for many application in the industry.

## How to get better results from the CNN
tweaking the differents parameters of the model could lead to better performance, changing the learning rate to find better weight that could help achieving a better accuracy, also it can be observe while looking at the training and validation error that the model is maybe is not overfittin at all and therefore might benefit of adding one additionnal layer, or changing the non linear function used in each layer. Overall, it is possible to achieve better results by tweaking the parameter and trial and error

increasing the number of data in the training dataset,  

using real measurement ??

using an approximate of the frequency response of a microphone ? 
## Room simulation parameters that can help achieving a better results 

As observed in the predicted DOA from the CNN, the model would probably perform better by increasing the angle resolution, however it would also mean using much more training exemple to generalize better in any direction.

Using a larger range of T60 to for the model to be applicable in more cases (when reverb is small or large) 

to improve the accuracy of the model on some specific data, it would be possible to simulate rooms with lower reverberation times as it is known that the reverberation in the measurement leads to more error. 
## Other dataset that could be used in order to train the ANN
In order to predict the direction of arrival it would we could imagine using a dataset of multichannel room impulse response.

