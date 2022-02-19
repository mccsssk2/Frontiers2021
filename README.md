# pm3 Lumped Parameter Human Circulation Code 2021
License:
GNU GENERAL PUBLIC LICENSE Version 3

2nd February 2021.

Please cite this reference: Joseph, J.J.; Hunter, T.J.; Sun, C.; Goldman, D.; Kharche, S.R.; McIntyre, C.W. Using a Human Circulation Mathematical Model to Simulate the Effects of Hemodialysis and Therapeutic Hypothermia. Appl. Sci. 2022, 12, 307. https://doi.org/10.3390/app12010307

DOI: ( https://doi.org/10.3390/app12010307).

Code description:
This code contains lumped parameter model capable of simulating whole body blood flow with detailed kidney blood flow, dialysis, and baroreflex control mechanism. The model
works as a system of stiff ODEs that is solved using MATLAB (2019b onwards) solver ode15s.

Inputs: The model itself is dependent on the parameters in p() which contains the many physical constants and physiological parameters of the model, as well as the initial
conditions of the state variables in y0(). These are not standard inputs, as they are hard coded, but can be modified to fit the needs of the user. Other capabilities
include modifying the parameters over time using the modParam function. This feature is activated by setting control = false, currently the modParam function activates a
sodium profile for the dialyzer unit.

Outputs: This code prints all outputs to a data file named output*.dat with * as the integer n. The output is arranged as follows: column 1 is time, columns 2-45 are the
state variables y() in order, columns 46 to 192 are the model parameters p() in order starting from p(2) (as p(1) is time).

References.

Heldt, T., Shim, E. B., Kamm, R. D., Mark, R. G., & Massachusetts. (2002). Computational modeling of cardiovascular response to orthostatic stress. Journal of
Applied Physiology, 92(3), 1239–1254. https://doi.org/10.1152/japplphysiol.00241.2001

Ursino, M., Colí, L., Brighenti, C., Chiari, L., De Pascalis, A., & Avanzolini, G. (2000). Prediction of solute kinetics, acid-base status, and blood volume
changes during profiled hemodialysis. Annals of Biomedical Engineering, 28(2), 204–216. https://doi.org/10.1114/1.245 

Lin, J., Ngwompo, R. F., & Tilley, D. G. (2012). Development of a cardiopulmonary mathematical model incorporating a baro–chemoreceptor reflex control system.
Proceedings of the Institution of Mechanical Engineers, Part H: Journal of Engineering in Medicine, 226(10), 787-803. doi:10.1177/0954411912451823

Levy, M. N., & Zieske, H. (1969). Autonomic control of cardiac pacemaker activity and atrioventricular transmission. Journal of Applied Physiology, 27(4), 465-
470. doi:10.1152/jappl.1969.27.4.465
