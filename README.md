# SpasTract
  Welcome to this repository dedicated to the SpasTract research project tutorial along official documents regarding the use of animal models.
SpasTract stands for "Spas(tin) Tract(ography)", our team thus try to elucidate the consequences of SPASTIN gene alteration (i.e in Hereditary Spastic Paraplegia or HSP) and ways to characterize it first on murine animal models, ex vivo/in vivo, and then later on we could shift discoveries to the human and allow for potential early diagnosis of HSP for example. 
Our team is looking at this topic through the scope of Magnetic Resonance Imaging or MRI, more particularly dMRI (diffusion MRI), which allows us in turn to do Tractography, the reconstruction of white matter fibers non-invasively. Reconstructed fibers are called 'streamlines' as they do not represent the ground truth but rather an indirect metric of the actual amount of fibers in the brain.
Our team also makes use of other tools to further complete this framework such as Histological cryocuts with gold marking, Carbocyanine dye DiI for axonal tracing, Compressed Sensing for high resolution MRI with partial acquisition of the K-space,... On top of that we also do a kind of Tracto'metry' based on Tractography with for example structural connectivity matrices.

Below you will find examples of work done by our team and that can be reproduce by following this tutorial :

![SpasTractSPASTINvsWT_Panel](https://github.com/user-attachments/assets/53911853-04f9-476a-a028-8dd538cb58d0)
  
#Examples of some major tracts : anterior commissure and corpus callosum both from a WT sample and then from a KO-SPASTIN sample.



![P56AtlasAnnotationProperOrientation](https://github.com/user-attachments/assets/401fcace-f0e5-4d1d-a367-7e4486698164)  
#An example of brain atlas availabe online from AllenBrain.

![MatrixMakingfromConnectome_STR-MB](https://github.com/user-attachments/assets/60a2ae10-cc22-42ef-9988-2b408cbef0fb) 
#An example of 'Connectomes' or Connectivity Matrices (structural) showcasing the amount of streamlines between regions based on Tractography co-registered with the atlas shown above.

![Exemple_REOR_3DCSSPASTINoldF1023A](https://github.com/user-attachments/assets/a3f90582-fc8d-400c-b111-a3d92dd7d68e)  
#An example of Compressed Sensing reconstruction based on undersampled K-space MRI acquisition (Acceleration Factor of 2 which means 50µm resolution can be done for the time of a 100µm resolution scan).

![VolumetrySignificantGMWM_V2](https://github.com/user-attachments/assets/437f38e2-3ab1-447a-b4cc-c1ca8f896c7b)
#An example of volumetry done based on an atlas coregistered to a high resolution Compressed Sensing image.

------------------------------
#### SpasTract Tutorial 
  You will find in the associated SpasTractTutorialFolder folder the SpasTract tutorial elements as well as a detailed step-by-step guide... 

