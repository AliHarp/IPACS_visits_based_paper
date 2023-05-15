# Planning integrated care services using simulation modelling for visits-based home care with time varying resource requirements

Alison Harper, David Worthington, Zehra Önen-Dumlu,  Paul Forte, Christos Vasilakis, Martin Pitt Richard Wood


PenCHORD, University of Exeter Medical School Rd Exeter EX1 1TE, UK
a.l.harper@exeter.ac.uk 	m.pitt@exeter.ac.uk

Department of Management Science, Lancaster University, Bailrigg, Lancaster LA1 4YW, UK
d.worthington@lancaster.ac.uk 

University of Bath School of Management Claverton Down Bath BA2 7AY, UK
zod23@bath.ac.uk c.vasilakis@bath.ac.uk 

NHS Bristol, North Somerset and South Gloucestershire CCG South Plaza, Marlborough Street Bristol BS1 3NX, UK
richard.wood16@nhs.net paul.forte@nhs.net 

Correspondence to: Alison Harper a.l.harper@exeter.ac.uk 

**Abstract**

Delays to hospital discharge impact the wider healthcare system and are a contributory factor to emergency department overcrowding, ambulance offload and response delays experienced both in the UK and globally. Health and care services comprising of care provided in community facilities and at home support timely discharge from hospital for patients with complex care needs. It attempts to avoid costly discharge delays and consequent risks to patient physical and mental health. Home care is provided by visiting care staff in patients’ homes, with time-varying visit requirements, typically reducing over time. This paper reports the development of a novel discrete-time stochastic simulation model in a large health and social care system in England, which models patients using a time-varying resource draw, rather than a fixed capacitated resource. Our results are compared to steady-state and time-dependent analytical approximations to identify areas of the parameter space in our case study where the analytical models perform sufficiently well, and to enable an empirical case for the use of time-varying resources to model visits-based home care. The flexibility of our simulation model and our purposeful engagement with stakeholders have enabled the resulting tool to be in routine use for home care resource planning. It is available open source for take-up by researchers or analysts involved in home care demand and capacity planning.

## Author ORCIDs

[![ORCID: Harper](https://img.shields.io/badge/ORCID-0000--0001--5274--5037-brightgreen)](https://orcid.org/0000-0001-5274-5037)
[![ORCID: Worthington](https://img.shields.io/badge/ORCID-0000--0003--1400--0194-brightgreen)](https://orcid.org/0000-0003-1400-0194)
[![ORCID: OnenDumlu](https://img.shields.io/badge/ORCID-0000--0001--8878--5495-brightgreen)](https://orcid.org/0000-0001-8878-5495)
[![ORCID: Forte](https://img.shields.io/badge/ORCID-0000--0002--1060--9106-brightgreen)](https://orcid.org/0000-0002-1060-9106)
[![ORCID: Wood](https://img.shields.io/badge/ORCID-0000--0002--3476--395X-brightgreen)](https://orcid.org/0000-0002-3476-395X)
[![ORCID: Vasilaki](https://img.shields.io/badge/ORCID-0000--0002--0391--0910-brightgreen)](https://orcid.org/0000-0002-0391-0910)
[![ORCID: Pitt](https://img.shields.io/badge/ORCID-0000--0003--4026--8346-brightgreen)](https://orcid.org/0000-0003-4026-8346) 
 
## About the IPACS Project

The IPACS project was funded by Health Data Research UK. It ran from May 2020-March 2023. The project was undertaken by a research team of six from Bristol, North Somerset, South Gloucestershire Integrated Care Board (BNSSG ICB) (Dr Richard Wood and Dr Paul Forte), University of Bath School of Management (Prof. Christos Vasilakis and Dr Zehra Onen Dumlu) and University of Exeter Medical School (Prof. Martin Pitt and Dr Alison Harper).

The IPACS project aimed to investigate what might constitute 'optimal capacity' along different parts of the complex care discharge pathways form acute hospital to community healthcare. 

The IPACS simulation model is a high-level computer model of the discharge-to-assess (D2A) pathways.  

## Documentation

The model is reported using [STRESS-DES Reporting Guidelines](https://doi.org/10.1080/17477778.2018.1442155) (Monks et al. 2019)

## Results 

The models and input files used for this paper are available here:

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.7937061.svg)](https://doi.org/10.5281/zenodo.7937060)

All model outputs and results are available here:

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.7934721.svg)](https://doi.org/10.5281/zenodo.7934721)

## Using and citing the model
The IPACS model is licensed using [GPL-3 license](https://choosealicense.com/licenses/gpl-3.0/).

If you use the IPACS routine model for research, reporting, education or any other reason, please cite it using details on Zenodo
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.7845908.svg)](https://doi.org/10.5281/zenodo.7845908)

```
@software{alison_harper_2023_7845908,
  author       = {Alison Harper and
                  Zehra Onen Dumlu and
                  Paul Forte and
                  Christos Vasilakis and
                  Martin Pitt and
                  Richard Wood},
  title        = {{Code and material for IPACS Model - Improving the 
                   Flow of Patients between Acute, Community and
                   Social Care}},
  month        = apr,
  year         = 2023,
  publisher    = {Zenodo},
  version      = {v1.0.0},
  doi          = {10.5281/zenodo.7845908},
  url          = {https://doi.org/10.5281/zenodo.7845908}
}
```




## License  
[GPL-3 license](https://choosealicense.com/licenses/gpl-3.0/)

Permissions of this strong copyleft license are conditioned on making available complete source code of licensed works and modifications, which include larger works using a licensed work, under the same license. Copyright and license notices must be preserved. Contributors provide an express grant of patent rights. 
## Repo overview

```bash

├── MODEL_steady_state.R
├── MODEL_time_dependent.R
├── images
│   └── visit_flow.png
├── intelligent_init.R
├── save_outputs_to_single_file.R
├── STRESS_DES
├── LICENSE
├── README.md
└── .gitignore
└── visit_vector_test.R
└── time_dependent_inputs.xlsx
└── time_dependent_inputs_E.xlsx
└── time_dependent_inputs_N.xlsx
```


## Research Papers linked to IPACS

* [**The False Economy of Seeking to Eliminate Delayed Transfers of Care: Some Lessons from Queueing Theory**](https://link.springer.com/article/10.1007/s40258-022-00777-2)

This study aims to demonstrate how, counter to intuition, pursual of elimination of acute delayed transfers of care is likely to be uneconomical, as it would require large amounts of community capacity to accommodate even the rarest of demand peaks, leaving much capacity unused for much of the time.

 
*  [**Optimising the balance of acute and intermediate care capacity for the complex discharge pathway: Computer modelling study during COVID-19 recovery in England**]( https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0268837)
 	
 A simulation study using the IPACS model on the effects of COVID-19 on community capacity requirements and costs to minimise acute delayed discharges.
 	
 
*  [**A Demand and Capacity Model For Home-Based Intermediate Care: Optimizing The'Step Down’Pathway**](https://ieeexplore.ieee.org/abstract/document/9715468)

An application of the discrete-time simulation model showing that total costs across the acute-community interface can be minimized by identifying optimal community capacity in terms of the maximum number of patients for which home visits can be provided by the service.


* [**Improving Hospital Discharge Flow Through Scalable Use of
Discrete Time Simulation and Scenario Analysis**](https://doi.org/10.36819/SW23.013)

This paper reports on the development and deployment of the versatile IPACS simulation tool for modelling both the home-based and bedded community step-down pathways, known as ‘Discharge to Assess’  in England’s NHS. Developed in open source ‘R’, the tool offers scalable solutions for exploring different scenarios relating to demand, capacity and patient length of stay. 

Additional information is available here:
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.7845995.svg)](https://doi.org/10.5281/zenodo.7845995) 

