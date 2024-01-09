# Human-Like Behavior of Automated Driving Systems
This repository contains algorithms for driver model analysis and prototypes for ADAS functions endowed with human-like features.

## How is this repository organized?
The most up-to-date versions of exact functional solutions can be found in the 'Solutions' subfolder. This contains further sub-folders with major solutions in the chain of the automated driving system.
The code contained by these folders, or the links to other code spaces are validated, their validation details can be found in their corresponding read-me files.
Further code can be found in the matlab_evaluation and the python_evaluation subfolders. These contain a more complex development chain, the code inside is still under development.

## Published Papers
[1] Ignéczi G., Horvath E., Pup D.: Implementation of a self-developed model predictive control scheme for vehicle parking maneuvers. The first ISTRC Annual Conference, Tel-Aviv, Izrael, 21, June 2021. \
https://www.researchgate.net/publication/354774900_Implementation_of_a_self-developed_model_predictive_control_scheme_for_vehicle_parking_maneuvers \
[2] G. Igneczi and E. Horvath: A Clothoid-base Local Trajectory Planner with Extended Kalman Filter. Paper presented at the Applied Machine Intelligence and Informatics (SAMI 2022) conference, Poprád, Slovakia, 2-5 March 2022\
https://ieeexplore.ieee.org/document/9780857 \
[3] G. Igneczi: Lane Keep Assist Optimized to Human Like Behavior. - Digital Vehicle Industry Researches at the University of Gyor, Gyor, Hungary, 15-17 June, 2022. \
https://jkk-web.sze.hu/esemeny/digitalis-jarmuipari-kutatasok-2022-mesterseges-intelligencia-a-mobilitasban \
[4] G. Igneczi: Validation of a Driver Model in a Closed-loop Simulation Environment - Conference of Mobility and Environment, Gyor, Hungary, 2-4. November, 2022.\
https://jkk-web.sze.hu/jovoformalo-jarmuipari-kutatasok-konferencia-2022/ \
[5] G. Igneczi and E. Horvath: Node Point Optimization for Local Trajectory Planners based on Human Preferences. IEEE 21st World Symposium on Applied Machine Intelligence and Informatics, Herlany, Slovakia, 19-21. January, 2023. \
https://ieeexplore.ieee.org/document/10044488 \
[6] Igneczi, G., Horvath, E., Toth, R., & Nyilas, K. (2023). Curve Trajectory Model for Human Preferred Path Planning of Automated Vehicles. Automotive Innovations - *in press* \
https://arxiv.org/ftp/arxiv/papers/2310/2310.02696.pdf \
[7] G. Igneczi, T. Dobay: Typification of Driver Models Using Clustering Methods - Conference of Sustainable Vehicle Industry Technologies at University of Gyor, Gyor, Hungary, 14-16, June, 2023. \
https://jkk-web.sze.hu/esemeny/fenntarthato-jarmuipari-technologiak-kutatasa-a-szechenyi-istvan-egyetemen-konferencia/ \
[8] G. Igneczi, E. Horvath and K. Nyilas: A Linear Driver Model of Local Path Planning for Lane Driving, Paper presented at IEEE 21st International Symposium on Intelligent Systems and Informatics (SISY 2023), Pula, Croatia, 21-23, September, 2023. \
http://conf.uni-obuda.hu/sisy2023/ \
[9] G. Igneczi, D. Jozsa: The Implementation of a Human-Like Lane Keep System in a Prototype Vehicle, Conference of Sustainable Vehicle Industry Technologies at University of Gyor, Gyor, Hungary, 2-4, November, 2023. \
https://jkk-web.sze.hu/fenntarthato-jarmuipari-technologiak-kutatasa-a-szechenyi-istvan-egyetemen-konferencia-2023-november-06-07/ \
[10] G. Igneczi, E. Horvath: Parameter Identification of the Linear Driver Model, in review \
\
[11] G. Igneczi, E. Horvath, A. Borsos and K. Nyilas: Short-term Compensation Strategy of Human Drivers in Lane Following, in review \
\

## Main Sources
### Literature of path planning (PPxx)
[PP1]   S. Alqahtani, I. Riley, S. Taylor and R. Mailler, "Predictive Path Planning Algorithm Using Kalman Filters and MTL Robustness," in IEEE International Symposium on Safety, Security, and Rescue Robotics, University of Pennsylvania, 2018.\
[PP2]   P. D. Ganesha, B. Subathra, G. Saravanakumar and S. Seshadhri, "Extended Kalman Filter Based Path-Planning Algorithm for Autonomous Vehicles with I2V Communication," in 4th IFAC Conference on Advances in Control and Optimization of Dynamical Systems ACODS 2016, Tiruchirappali, India, 2016, February. \
[PP3]   B. Amol, H. Monson and S. T. Mark, "Robust lane detection and tracking with ransac and Kalman filter," in 16th IEEE International Conference on Image Processing (ICIP), Cairo, Egypt, 2009. \
[PP4]   D. Byambaa, H. Sabir and L. Deok-Jin, "Highly Curved Lane Detection Algorithms," Applied Sciences, vol. X., 2020. \
[PP5]   D. H. Douglas and T. K. Peucker, "Algorithms for the reduction of the number of points required to represent a digitized line or its caricature," Cartographica: the international journal for geographic information and geovisualization, vol. 10, no. 2, pp. 112-122, 1973. \
[PP6]   E. Bertolazzi and M. Frego, "G1 Fitting with Clothoids," Wiley Online Library, 7 March 2014. \
[PP7]   C. Gackstatter, T. Sven , P. Dr. Heinemann and G. Prof. Klinker, "Stable Road Lane Model Based on Clothoids," in Advanced Microsystems for Automotive Applications, Springer, 2010, pp. 133-143. \
[PP8]   J. Xiao, L. Luo, Y. Yao, W. Zou and R. Klette, "Lane Detection Based on Road Module and Extended Kalman Filter," 2018. \
[PP9]   P. Pokojska and W. Pokojski, "Voronoi diagrams – inventor, method, applications," Polish Cartographical Review, vol. 50, no. 3, pp. 141-150, 2018. \
[PP10]  A. Elfes, "Using Occupancy Grids for Mobile Robot Perception and Navigation," in COMPUTER, Pittsburgh, Carnegie Mellon University, 1989, pp. 46-57. \
[PP11]  D. Gonzalez-Arjona, A. Sanchez, F. Lopez-Colino, A. de Castro and J. Garrido, "Simplified Occupancy Grid Indoor Mapping Optimized for Low-Cost Robots," ISPRS International Journal of Geo-Information, vol. 2, pp. 959-977, 2013. \
[PP12]  Y. Jiamin and W. Hao, "Research on a Costmap that can Change the Turning Path of Mobile Robot," Journal of Physic, pp. 1-17, 2021. \
[PP13]  M. McNaughton, C. Urmson, J. M. Dolan and J.-W. Lee, "Motion planning for autonomous driving with a conformal spatiotemporal lattice," in IEEE International Conference on Robotics and Automation (ICRA) , Shanghai, China, 2011.\
[PP14]  N. Winston, "Continuous-Curvature Paths for Autonomous Vehicles," in Proceedings, 1989 International Conference on Robotics and Automation, Scottsdale, USA, 1989.\
[PP15]  Y. Wang, D. Shen and K. E. Teoh, "Lane detection using spline model," Pattern Recognition Letters, vol. 21, no. 8, pp. 677-689, 2000.\
[PP16]  C. R. Pozna and E. Horvath, "Clothoid-based Trajectory Following Approach for Self-driving vehicles," in IEEE 19th World Symposium on Applied Machine Intelligence and Informatics (SAMI), Herlany, Slovakia, 2021.\
[PP17]  M. Fatemi, L. Hammarstrand, L. Svensson and G.-F. Ángel F." Road Geometry Estimation Using a Precise Clothoid Road Model and Observations of Moving Vehicles," in 17th International IEEE Conference on Intelligent Transportation Systems (ITSC), Qingdao, China , 2014.\
[PP18]  SAFESTAR, "Safety Standards for Road Design and Redesign," European Comission, 2002.\
[PP19]  B. Enrico and F. Marco, "Fast and accurate clothoid fitting," Optimal control applied to human motion and physiology, September, 2012.\

### Literature of traditional driver modelling (Dxx)
[D1] 	E. A. Fleishman, „Development of a behavior taxonomy for describing human tasks: A correlational-experimental approach,” Journal of Applied Psychology, pp. 1-10, 1967. \
[D2] 	J. A. McKnight és B. B. Adams, „Driver Education Task Analysis,” Human Resources Research Organisation, Alexandria, Va., 1970. \
[D3] 	L. Evans és R. C. Schwig, „Human behavior and traffic safety,” Plenum Press, pp. 485-520, 1985. \
[D4] 	G. J. S. Wilde, „The Theory of Risk Homeostasis: Implications for Safety and Health,” Risk Analysis, (1), vol. 4, pp. 209-225, 1982. \
[D5] 	D. Klebersberg, „ Subjektive und objektive Sicherheit im Strassenverkehr als Aufgabe für die Verkehrssicherheitsarbeit,” Schriftenreihe der Deutschen Verkehrswacht, (1), pp. 3-12, 1971. \
[D6] 	R. Fuller, „A conceptualization of driving behaviour as threat avoidance,” Ergonomics, (27), vol. 1, pp. 1139-1155, 1984. \
[D7] 	D. C. Wickens, „Engineering Psychology and Human Performance,” Foresman and Company, Glenview, Illinios, 1984. \
[D8] 	J. Theeuwes, „Visual attention and driving behaviour,” in In Proceedings of the International Seminar on Human Factors in Road Traffic, Braga, Portugal, 1993. \
[D9] 	C. Delphine és T. Gordon, „TRB Workshop on Driver Models: A Step Towards a Comprehensive Model of Driving,” Modelling Driver Behaviour in Automotive Environments, pp. 26-42, 2007. \
[D10] 	M. Panou, E. Bekiaris és V. Papakostopoulos, „Modelling Driver Behaviour in European Union and International Projects,” Modelling Driver Behaviour in Automotive Environments, pp. 3-25, 2007. \
[D11] 	K. Sarvesh, d. W. Joost és A. David, „Human-like driving behaviour emerges from a risk-based driver model,” Nature Communications, pp. 4850-4863, 2020. \
[D12]   Nicola Bongiorno, Orazio Pellegrino, et al. Cumulative lateral position: a new measure for driver performance in curves. Traffic Safety Research, 4(1):20–34, 2022. \
[D13]   Hongsheng Qi and Xianbiao Hu. Behavioral investigation of stochastic lateral wandering patterns in mixed traffic flow. Transportation Research Part C, 4(1):1–26, 2023. \
[D14]   Gao Hongbo, Xie Guotao, et al. Lateral control of autonomous vehicles based on learning driver behavior via cloud model. The Journal of China Universities of Posts and Telecommunications, 24(2):10–17, 2017.

### Literature of control driver models (Cxx)
[C1] 	F. Mars és C. Philippe, „Modelling human control of steering for the design of advanced driver assistance systems,” in Annual Reviews in Control, 2017. \
[C2] 	U. Kiencke, R. Majjad és S. Kramer, „Modeling and performance analysis of a hybrid driver model,” Control Engineering Practice, (7), vol. 8, pp. 985-991, 1999. \
[C3] 	C. Rathgeber, F. Winkler, D. Odenthal és S. Müller, „Lateral trajectory tracking control for autonomous vehicles,” in European Control Conference (ECC), Strasbourg, France, 2014. \
[C4] 	D. D. Salvucci és R. Gray, „A two-point visual control model of steering,” Perception, (1), vol. 33, pp. 1233-1248, 2004. \
[C5] 	C. R. Conlter, Implementation of the Pure Pursuit Path Tracking Algorithm, Pittsburgh, Pennsylvania: The Robotics Institute Carnegie Mellon University, 1992. \
[C6] 	C. C. McAdam, „An Optimal Preview Control for Linear Systems,” Journal of Dynamic Systems, Measurement and Control, (102), vol. 1, pp. 188-190, 1980. \
[C7] 	A. Y. Ungoren és H. Peng, „An Adaptive Lateral Preview Driver Model,” Vehicle System Dynamics, (1), vol. 43, pp. 245-259, 2005. \
[C8] 	P. Björn és L. Nilsson, „Modelling the Driver in Control,” Modelling Driver Behaviour in Automotive Environments, pp. 85-104, 2007. \
[C9] 	T. Choudhari és A. Maji, „Modeling driver’s braking and steering behavior along horizontal curves of two-lane rural highways for ADAS applications,” Traffic Injury Prevention, (7), vol. 23, pp. 404-409, 2022. \
[C10]   Horvath Erno, Hajdu Csaba, and Koros Peter. Novel pure-pursuit trajectory following approaches and their practical applications. In Proceedings of the 10th IEEE International Conference on Cognitive Infocommunications, pages 1–6, Naples, Italy, 2019. \
[C11]   R.A. Hess and A. Modjtahedzadeh. A control theoretic model of driver steering behavior. IEEE Control Systems Magazine, 10(5):3–8, 1990. \
[C12]   Haobin Jiang, Huan Tian, and Yiding Hua. Model predictive driver model considering the steering characteristics of the skilled drivers. Advances in Mechanical Engineering, 11(3):1–14, 2019. \
[C13]   Alexander Katriniok, Jan P. Maschuw, et al. Optimal vehicle dynamics control for combined longitudinal and lateral autonomous vehicle guidance. In Proceedings of European Control Conference, pages 974–979, Z¨urich, Switzerland, 2013. \
[C14]   J. Kong, M. Pfeiffer, G. Schildbach and F. Borelli, "Kinematic and Dynamic Vehicle Models for Autonomous Driving Control Design," in IEEE Intelligent Vehicles Symposium, Seoul, Korea, 2015.\
[C15]   H. Ye, H. Jiang, S. Ma, B. Tang and L. Wahab, "Linear model predictive control of automatic parking path tracking with soft constraints," International Journal of Advanced Robotic Systems, pp. 1 - 13, May - June 2019. \
[C16]   O. Pauca, C. F. Curuntu and C. Lazar, "Predictive Control for the lateral and longitudinal dynamics in automated vehicles," in 23rd International Conference on System Theory, Control and Computing, Sinaia, 2019. \
[C17]   L. Wang, Model Predictive Control System Design and Implementation using Matlab, London: Springer - Verlag, 2009, pp. 28 - 84, 324 - 359.

### Literature of planner driver models (Pxx)
[P1] 	I. Bae, J. Moon és J. Seo, „Toward a Comfortable Driving Experience for a Self-Driving Shuttle Bus,” Electronics, (8), vol. 9, pp. 943-955, 2019. \
[P2] 	K. Christos, Q. Mohammed, C. Wen-Hua és D. Lipika, „Real-time motion planning methods for autonomous on-road driving: State-of-the-art and future research directions,” Transportation Research Part C: Emerging Technologies,(1), vol. 60, pp. 416-442, 2015. \
[P3] 	F. Braghin, F. Cheli, S. Melzi és E. Sabbioni, „Race Driver Model,” Computer and Structures, %1. szám86, pp. 1503-1516, 2008. \
[P4] 	Y. Cenxin, N. Anning, L. Jing , W. Jinghui, Z. Chunqin, C. Qinqin és T. Yifeng, „A Novel Dynamic Lane-Changing Trajectory Planning Model for Automated Vehicles Based on Reinforcement Learning,” Journal of Advanced Transportation, (1), vol. 1, pp. 1-16, 2022. \
[P5] 	T. Gu és J. M. Dolan, „Toward human-like motion planning in urban environments,” in IEEE Symposium on Intelligent Vehicle, Dearborn, USA, 2014. \
[P6] 	P. Iakovos és T. Masayoshi , „Fast Lane Changing Computations using Polynomials,” in Proceedings of the American Control Conference, Denver, USA, 2003. \
[P7] 	E. Lambert, R. Romano és D. Watling, „Optimal Path Planning with Clothoid Curves for Passenger Comfort,” in 5th International Conference on Vehicle Technology and Intelligent Transport Systems, Heraklion, Crete, 2019. \
[P8] 	W. Moritz, Z. Julius, K. Sören és T. Sebastian, „Optimal trajectory generation for dynamic street scenarios in a Frenét Frame,” in IEEE International Conference on Robotics and Automation, 2010. \
[P9] 	O. Takayuki, „Multimodal trajectory optimization for motion planning,” The International Journal of Robotics Research,(39), vol. 8, pp. 983-1001, 2020. \
[P10] 	O. Siebinga, A. Zgonnikov és A. David, „A Human Factors Approach to Validating Driver Models for Interaction-aware Automated Vehicles,” ACM Transactions on Human-Robot Interaction ,(11), vol. 4, pp. 1-21, 2022. \
[P11] 	L. Xiao, L. Jun és Z. Hua, „Dynamic motion planner with trajectory optimisation for automated highway lane-changing driving,” IET Intelligent Transport Systems, (14), vol. 14, pp. 2133-2140, 2021. \
[P12] 	H. Xu, X. Donghao, Z. Huijing, M. Mathieu, A. Francois és G. Frank, „A Human-like Trajectory Planning Method by Learning from Naturalistic Driving Data,” in IEEE Intelligent Vehicles Symposium (IV), Suzhou, China, 2018. \
[P13] 	Kazuhide, O. Kazuhide, I. Laurent és T. Panagiotis, „Vision-Based Autonomous Path Following Using a Human Driver Control Model With Reliable Input-Feature Value Estimation,” IEEE Transactions on Intelligent Vehicles, pp. 497-506, 2019. \

### Literature of human and ADAS (HUxx)
[HU1] 	M. M. Aydin, „A new evaluation method to quantify drivers’ lane keeping behaviors on urban roads,” Transportation Letters, (10), vol. 12, pp. 738-749, 2020. \
[HU2] 	Z. Yassine, E. Khalid, R. Salahheddine, B. Afaf és B. Ismail, „Driver profiling: The pathway to deeper personalization,” Journal of King Saud University - Computer and Information Sciences, (10), vol. 34, pp. 9088-9101, 2022. \
[HU3] 	M. Zanker, F. Ricci, D. Jannach és L. Terveen, „Measuring the impact of personalization and recommendation on user behaviour,” International Journal of Human-Computer Studies, (8), vol. 68, pp. 469-471, 2010. \
[HU4] 	R. E. Chandler, R. Herman és E. W. Montroll, „Traffic Dynamics: Studies in Car Following,” Operations Research, pp. 165-184, 1958. \
[HU5] 	W. Zheng, K. Tsutomu és N. Kimihiko, „Effect of Haptic Guidance Steering on Lane Following Performance by Taking Account of Driver Reliance on the Assistance System,” in 2018 IEEE International Conference on Systems, Man, and Cybernetics (SMC), Miyazaki, Japan, 2018. \
[HU6] 	A.-C. Hensch, M. Beggiate és J. F. Krems, „Drivers’ gap acceptance during parking maneuvers as a basis for initiating driving actions in automated vehicles,” Transportation Research Part F: Traffic Psychology and Behaviour, pp. 133-142, 2023. \
[HU7] 	A.-C. Hensch, S. Mandl, M. Beggiato és S. Anja, „The interplay of personality traits with drivers’ gap acceptance,” in 13th International Conference on Applied Human Factors and Ergonomics (AHFE 2022), New York, USA, 2022. \
[HU8] 	Z. K. Zhaobo, K. Akash és M. Teruhisa, „Detection of Perceived Discomfort in SAE L2 Automated Vehicles through Driver Takeovers and Physiological Spikes,” in 2022 IEEE 25th International Conference on Intelligent Transportation Systems (ITSC), Macau, China, 2022.

### Definition of ADAS (Axx)
[A1] 	SAE International, „Adaptive Cruise Control (ACC) Operating Characteristics and User Interface,” 2014. \
[A2] 	B. Nguyen, N. Famiglietti, O. Khan és R. Hoang, „Testing and Analysis of Lane Departure Warning and Lane Keeping Assist System Response,” SAE Int. J. Adv. & Curr. Prac. in Mobility, (3), vol. 5, pp. 2301-2316, 2021. \
[A3] 	SAE International, „Taxonomy and Definitions for Terms Related to Driving Automation Systems for On-Road Motor Vehicles,” 2021. \
[A4] 	D. G. Bautista, Function architecture for automated vehicles trajectory planning in complex environments, Paris: Université Paris sciences et lettres, 2017. \
[A5] 	R. Novickis, R. Kadikis, A. Levinskis és V. Fescenko, „Functional Architecture for Autonomous Driving and its Implementation,” in Conference of 17th Biennial Baltic Electronics Conference (BEC), Tallin, Estonia, 2020.

### Clustering (CLxx)
[CL1] 	J. Hartigan, Clustering algorithms, New York, NYUnited States: John Wiley & Sons, 1975. \
[CL2] 	K. Sasirekha és P. Baby, „Agglomerative Hierarchical Clustering Algorithm - A Review,” International Journal of Scientific and Research Publications, (3), vol. 3, pp. 1-3, 2013. \
[CL3] 	M. Cui, „Introduction to the K-Means Clustering Algorithm Based on the Elbow Method,” Geoscience and Remote Sensing, (3), vol. 3, pp. 9-16, 2020. \
[CL4] 	L. Kaufman és P. J. Rousseeuw, Finding Groups in Data: An Introduction to Cluster Analysis, Hoboken, NJ: John Wiley & Sons, 1990. \
[CL5] 	J. S. Farris, „On the Cophenetic Correlation Coefficient,” Systematic Biology, (18), vol. 3, pp. 279-285, 1969.

### State Estimation (Exx)
[E1]    B. S´aez Pablo and Bruce E. Rittmann. Model-parameter estimation using least squares. Water Research, 26(6):789–196, 1992. \
[E2]    Da-Zheng Feng, Zhen Bao, and Li-Cheng Jiao. Total least mean squares algorithm. IEEE Transactions on signal processing, 46(8):1–10, 1998. \
[E3]    Radek Martinek, Jaroslav Rzidky, Rene Jaros, et al. Least mean squares and recursive least squares algorithms for total harmonic distortion reduction using shunt active power filter control. Energies, 12(8):1–10, 2019. \
[E4]    Paulo Sergio Diniz Ramirez. The least-mean-square (lms) algorithm. Part of the The Kluwer International Series in Engineering and Computer Science book series: Adaptive Filtering, 694(1):79–183, 1998. \
[E5]    Jorma Rissanen, Teemu Roos, and Petri Myllym¨aki. Model selection by sequentially normalized least squares. Journal of Multivariate Analysis, 101:839–849, 2010. \
[E6]    Neil J. Bershad and Jose C.M. Bermudez. Analysis of the least mean square algorithm with processing delays in the adaptive arm for gaussian inputs for system identification. Signal Processing, 11(8):1–10, 2023. \
[E7]    Ubaid M. Al-Saggaf, Muhammad Moinuddin, Muhammad Arif, et al. The q-least mean squares algorithm. Signal Processing, 11(1):50–60, 2015.\
[E8]    D. Semino, M. Moretta, and C. Scali. Parameter estimation in extended kalman filters for quality control in polymerization reactors. Computers and Chemical Engineering, 20(2):913–918, 1996.\
[E9]    Xiaodian Sun, Li Jin, and Momiao Xiong. Extended kalman filter for estimation of parameters in nonlinear state-space models of biochemical networks. PLoS ONE, 3(11):1–2, 2023.\
[E10]   Wooyoung Na and Chulsang Yoo. Real-time parameter estimation of a dual-pol radar rain rate estimator using the extended kalman filter. Remote Sens., 13(12):1–2, 2021.\
[E11]   Rudolph Emil Kalman. A new approach to linear filtering and prediction problems. Transactions of the ASME–Journal of Basic Engineering, 82(Series D):35–45, 1960.\
[E12]   Bruno Johannes Schnekenburger. A modified extended kalman filter as a parameter estimator for linear discrete-time systems. New Jersey Institute of Technology, Electronic Thesis and Dissertation, 1988.\
[E13]   Devyani Varshney, Mani Bhushan, and Sachin C. Patwardhan. State and parameter estimation using extended kitanidis kalman filter. Journal of Process Control, 76(1):98–111, 2019.\
[E14]   Yishen Zhao, Philippe Chevrel, et al. Continuous identification of driver model parameters via the unscented kalman filter. IFAC-Papers Online, 52(28):126–133, 2019.\
[E15]   Christina Gackstatter, Patrick Heinemann, et al. Stable road lane model based on clothoids. Conference Proceedings of Advanced Microsystems for Automotive Applications, pages 133–143, 2010.\
[E16]   Xin Li, Anzhi Lei, et al. Improving kalman filter for cyber physical systems subject to replay attacks: An attack-detection-based compensation strategy. Applied Mathematics and Computation, 466(128444):1–14, 2023.\
[E17]   R. Kleinbauer, "Kalman Filtering Implementation with Matlab," University Stuttgart Institute of Geodasy, Helsinki, 2004. \
[E18]   M. M. Shyam, N. Naren, R. O. Gemson and M. Ananthasayanam, "Introduction to the Kalman Filter and Tuning its Statistics for Near Optimal Estimates and Cramer Rao Bound," Department of Electrical Engineering Indian Institute of Technology, Kanpur, 2015. \
[E19]   M. K. Raman, "On the Identification of Variances and Adaptive Kalman Filtering," in IEEE Transactions on Automatic Control, 1970. \
[E20]   C. Cuadras, "On the Covariance between Functions," Journal of Multivariate Analysis, pp. 19-27, 2002.


# About the Author

Gergo Ferenc Igneczi has graduated as Control Engineer at Budapest University of Technology and Economics and learnt master studies in Vehicle Engineering at Széchenyi Istvan University of Győr, Hungary. He is a PhD student since 2021. As the part of his research, he studies the motion planning and control strategies of human drivers. Also, he is an active software engineer at Robert Bosch Kft, Budapest, Hungary since 2016, developing driver models and human-like motion planners and controllers for Advanced Driver Assistance Systems.

![alt text](profile.jpg)
