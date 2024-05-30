# 497R - Airfoil Analysis - Samuel Nasman

*Note: There is a wiki of unfamiliar terms in my 497R GitHub repository [here](https://github.com/samiamandyesican/497R/wiki).*

## Introduction and Methods

My goal in this project is to understand how an airfoil works, and how it is influenced under various parameters. For this project we used XFoil to explore the effects of angle of attack, thickness, and camber as well as Reynold's number. Additionally, we compared results to experimental data.

## Methods

The analysis of various parameters was airfoil performance was computed using a version of XFoil by Mark Drelas, adapted for use via the Julia programming language. Airfoil geometries were obtained from [airfoiltools.com](http://airfoiltools.com/airfoil/naca4digit). For exploring angle of attack and Reynold's number, the geometry of a NACA 2412 airfoil was used. Additionally, results from this geometry were compared to [experimental data from NASA](https://ntrs.nasa.gov/api/citations/19950002355/downloads/19950002355.pdf). To explore thickness and camber, results were generated for various NACA 4-series airfoils. Generally, analysis was done by calculating lift coefficient, drag coefficient, or moment coefficient over a range of angles of attack and then comparing results.  

## Results and Discussion
### Exploring Angle of Attack
One of the most obvious things to explore to understand how an airfoil works is to change angle of attack. Using the coordinates for a NACA 2412 airfoil with Reynold's number $2.2 \cdot 10^6$ and Mach number $0.13$, we obtain the following results from XFoil over a range of $-9^o \leq \alpha \leq 25^o $.

<img src="naca2412coor.png" alt="NACA 2412" width="49%" /> 
<img src="naca2412cl.png" alt="NACA 2412 Lift Coefficient" width="49%" /> 
<img src="naca2412cd.png" alt="NACA 2412 Drag Coefficient" width="49%" /> 
<img src="naca2412cm.png" alt="NACA 2412 Moment Coefficient" width="49%" /> 
<br><br>

***Lift Coefficient***  
As angle of attack increases from -5&deg; to around 12&deg;, we see a positive linear coorelation between angle of attack and lift coefficient. This is as expected since an increase in angle of attack leads to  
1) a decrease in the velocity of the air beneath the airfoil as it collides with the airfoil's lower surface (resulting in higher pressure) and 
2) an increase in the velocity above the airfoil as it is accelerated to fill the lower pressure wake of air displaced by the airfoil. 

This relationship that an increase in velocity is associated with a decrease in pressure is expressed by Bernoulli's principle. Bernoulli's equation makes several assumptions, but for our purposes the general coorelation between velocity and pressure still holds. 

$$
\red{P} + \frac{\rho}{2}\red{V^2}+\rho g h = constant
$$

At around 12&deg;, the air is too inertial to accelerate around the airfoil fast enough to fill the low pressure wake. This results in flow separation, also known as stall, (see image below) and rather than air being accelerated into the wake, a turbulent vacuum remains that increasingly pulls the airfoil backwards (increasing drag) rather than provide lift as the angle of attack continues to increase. 

![alt text](image-20.png)


***Drag Coefficient***  
As angle of attack increases there is an increase in the drag coefficient. This can be attributed to because there is a greater area of the lower surface of the airfoil normal to the freestream, and therefore more air pushing the airfoil backwards. The sharp increase in drag at high angles of attack can also be attributed to flow separation (stall) as described in the previous section.

***Moment Coefficient***  
We can see that this airfoil has a negative moment across angles of attack. This is common to most airfoils. Note that a negative moment for an airfoil means that it tends to pitch the aircraft down. This may be surpising since it seems intuitive that since the center of lift is normally at the quarter-chord, the airfoil would have a positive moment. However, if you look at the pressure distribution across an airfoil's surface (see image below), you can see that there is still negative pressure coefficient farther towards the trailing edge along the upper edge of the airfoil. This pressure has more influence on moment since $M= r \times F$ and there is a longer moment arm. 

![alt text](image-22.png)

As the airfoil approaches stall, the pressure distribution on the lower surface stays pretty much the same while on the upper surface as flow begins to separate the negative pressure becomes more skewed towards the leading edge as the separation point travels from the trailing edge towards the leading edge. This means the overall moment becomes less and less negative, approaching a maximum close to zero at stall and then quickly plunging negative as the pressure distribution across the entire upper surface becomes the same (see plots below. different airfoil, but same principle). In the plots below you can see that between alpha = 20 and alpha =24 the moment would become much more negative (as we can see in our plot as alpha exceeds 18 degrees or so)

![alt text](image-23.png)

So, generally we see that increasing angle of attack increases both lift and drag, up to the point of stall at which drag takes over. Therefore it's clear to see that an airfoil has a practical operating range or angles, outside of which the airfoil becomes pretty much useless (at least if your goal is to generate lift).

### Comparison to experimental data
Now I'm going to explore how well XFoil compares to real life. In other words I want to see the limits of XFoil so I can tell when it's useful and when it's inaccurate. I will first compare to NASA experimental data cited at the end of this document, using a Reynold's number 2.2e6 and Mach number 0.13. These two parameters are very important to ensure that the fluid is going to behave accurately in the simulation compared to the experiment because. Particularly Reynold's number since it describes the ratio between kinetic and viscous forces, which is the main factor in determining at what angle stall will occur. 

<img src="NASAvXFoil_cl.png" alt="NACA 2412 Lift Coefficient (XFoil v. NASA)" width="49%" /> 
<img src="NASAvXFoil_cd.png" alt="NACA 2412 Drag Coefficient (XFoil v. NASA)" width="49%" /> 
<img src="NASAvXFoil_cm.png" alt="NACA 2412 Moment Coefficient (XFoil v. NASA)" width="99%" /> 

From these images we can see that XFoil matches pretty well for lift coefficient and drag coefficient, although accuracy decreases at extreme angles of attack. This is unsurprising since the assumptions of panel theory break down when flow separates (the Kutta condition and irrotational flow for example, become unreasonable assumptions). The graph for moment coefficient v. lift coefficient however kind of seems all over the place. In the NASA paper it says that there were some discrepancies with their moment measurements and they weren't sure how to reconcile them, so perhaps it's not a problem with XFoil.  

Other experimental data by Dr. Abbot (sited below) uses a Reynolds number of 5.7e6. The graphs are below.
<img src="image-15.png" alt="ReSweepcm.png" width="49%" />
<img src="image-16.png" alt="ReSweepcm.png" width="49%" />
<img src="image-18.png" alt="ReSweepcm.png" width="49%" />

We can see once again that XFoil does a pretty good job but breaks down at high angles of attack. ADditionally we see that the moment diagram seems to match a little better, which is great! 

Overall we see that generally XFoil works best if we limit ourselves to reasonable angles of attack.


### Exploring Reynold's Number

Next I wanted to see how well XFoil handles a range of Reynolds' numbers. From the graphs you can see that there is siginificant numerical instability for very low Reynold's numbers. This makes sense since one of the assumptions we make in panel method is inviscid flow, which is associated with an infinitely high reynold's number. Based on the grphas I would be cautious going below 1e6 and so limited my graphs to that.


<img src="ReSweep.png" alt="ReSweep.png" width="49%" />
<img src="ReSweep2.png" alt="ReSweep2.png" width="49%" />
<img src="ReSweepcd.png" alt="ReSweepcd.png" width="49%" />
<img src="ReSweepcm.png" alt="ReSweepcm.png" width="49%" />

We can see pretty clearly from the graphs that increasing reynolds' number delays stall / flow separation (since a less inertial fluid is more easily able to accelerate around the airfoil without as much resistance). Therefore for high reynold's numbers, airfoils can use a wider range of angles of attack without stalling out. We also see that the airfoils have about the same slope before stall, which is a consequence of non-dimensionalization for the various coefficients. 

### Exploring Maximum Thickness
Next I wanted to see what happens if we change of the shape of the airfoil itself. How can it be optimized? First I looked at changing the thickness by generating different NACA 4-series (see graphs below).

<img src="image-1.png" alt="ReSweepcm.png" width="49%" />
<img src="image.png" alt="ReSweepcm.png" width="49%" />
<img src="image-13.png" alt="ReSweepcm.png" width="49%" />
<img src="image-2.png" alt="ReSweepcm.png" width="49%" />

As it turns out, XFoil really doesn't like airfoils that are extremely thin. This is likely because XFoil is making calculations very close to the center of the vortices used in panel method to simulate the shape of the airfoil. They have the velocity
$$
-\frac{\Gamma \times \hat{r}}{2 \pi r}
$$
and so as a result since computers aren't 100% precise there's going to be significant numerical instability as r approaches 0. 

Other than that thicker airfoils decrease lift and increase drag, which seems unfavorable. However they also have less abrupt stalls. Near the linear region all of the airfoils behave pretty much the same other than the largest airfoil. I would suspect this is because the airfoil is so thick that it causes some flow separation evean at low angles of attack. In general, increasing thickness decreases lift to drag ratio but also softens stall.


### Exploring Maximum Camber
Next I wanted to explore how changing the curviness (camber) changed airfoil performance. 


<img src="image-5.png" alt="ReSweepcm.png" width="49%" />
<img src="image-3.png" alt="ReSweepcm.png" width="49%" />
<img src="image-14.png" alt="ReSweepcm.png" width="49%" />
<img src="image-4.png" alt="ReSweepcm.png" width="49%" />

From the graphs it seems like increasing camber increases lift. However it also leads to more drag and eventually stall at lower angles of attack. This makes sense since as the airfoil is so curved, rotating it would make the lower surface more perpendicular to the free stream than if it weren't curved. Additionally the curved upper edge would be more likely to have flow separation at lower angles of attack. Looking at cl/cd we can see that for very low angles of attack, high camber generates a greater ratio of lift to drag, but for higher angles of attack lower camber tends to do better up to a certain point (around 20 degrees) at which they all pretty much converge. At this point both the least and most cambered have stalled out competely and so the shape doesn't do much difference. 

In general we can conclude that highly cambered airfoils increase lift, with drag and stall taking over at lower angles compared to a non-cambered airfoil. 

### Exploring Maximum Camber Position

Lastly I wanted to explore the effects of camber position on airfoil performance. I wasn't relly sure what to expect starting out, but based on the results its seems like it doesn't make much of a difference except for the moment. It turns out that a camber position farther towards the trailing edge increases the strength of the moment. This makes sense to me since the little "hook" at the end is so far from quarter chord, its flicking the air doward would have a larger effect since $M = r \times F$.

<img src="image-6.png" alt="ReSweepcm.png" width="49%" />
<img src="image-7.png" alt="ReSweepcm.png" width="49%" />
<img src="image-8.png" alt="ReSweepcm.png" width="49%" />
<img src="image-9.png" alt="ReSweepcm.png" width="49%" />
<img src="image-11.png" alt="ReSweepcm.png" width="49%" />
<img src="image-12.png" alt="ReSweepcm.png" width="49%" />

From these graphs, since camber position doesn't seem to make a huge difference on anything else, I'd want to minimize moment which still keeping it negative for consistency. Therefore having a camber position around a quarter chord seems most conventional.

## Conclusion

From the analyses and comparisons I've made I've come to a few conclusions about airfoil performance and the usefulness of XFoil. First, increasing angle of attack increases both lift and drag until stall, at which point drag takes over. Second, XFoil doesn't work very well for high angles of attack and low reynold's numbers. Third, a thinner airfoil generates more lift but also has more aggressive stall. Fourth, High camber generally creates more lift and is good for low angles of attack, but at high angles of attack there's a significant trade-off for higher drag. Fifth, camber position doesn't make a huge difference for anything except the moment, and upon analysis makes most sense to have around a quarter chord. 

## Works Cited
- http://www.airfoiltools.com/
- (nasa experimental data) https://ntrs.nasa.gov/api/citations/19950002355/downloads/19950002355.pdf
- webplot digitizer
- Abbott, I.H. and Von Doenhoff, A.E. (1959) Theory of Wing Sections: Including a Summary of Airfoil Data. Dover Publications, Mineola. https://www3.nd.edu/~ame40431/AME20211_2021/Other/AbbottDoenhoff_TheoryOfWingSectionsIncludingASummaryOfAirfoilData.pdf