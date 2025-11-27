
# Table of Contents

1.  [Project Description](#org3400e64)
2.  [Installation Instructions](#orgec3aa1b)
3.  [Usage Guide](#org4ff40c3)
4.  [Parameter Explanations](#org956e252)
5.  [Example Outputs](#org998d5b4)


<a id="org3400e64"></a>

# Project Description

"Simulation" is a project that measures the effect of behavior changes on the spread of COVID-19. It answers questions such as how going out less, socially distancing from others, self-isolating when sick, and reducing contact with others in general affects the spread of the COVID-19. The implementation consists of a stochastic SIR based model with behavioral component parameters.


<a id="orgec3aa1b"></a>

# Installation Instructions

To install dependencies for this program do the following:

$ python -m venv .venv
$ source .venv/bin/activate
$ pip install -r requirements


<a id="org4ff40c3"></a>

# Usage Guide

python simulation.py [-h] [-seed SEED] [-N<sub>S0</sub> N<sub>S0</sub>] [-N<sub>P0</sub> N<sub>P0</sub>] [-N<sub>A0</sub> N<sub>A0</sub>] [-N<sub>Y0</sub> N<sub>Y0</sub>] [-N<sub>R0</sub> N<sub>R0</sub>] [-s S] [-p P] [-a A] [-y Y] [-cycles CYCLES] [-avg<sub>steps</sub> AVG<sub>STEPS</sub>] plot<sub>name</sub> metrics<sub>name</sub> N<sub>S</sub><sub>name</sub> N<sub>P</sub><sub>name</sub> N<sub>A</sub><sub>name</sub> N<sub>Y</sub><sub>name</sub> N<sub>R</sub><sub>name</sub> time<sub>name</sub>


<a id="org956e252"></a>

# Parameter Explanations

positional arguments:
  plot<sub>name</sub>             plot file name to output
  metrics<sub>name</sub>          metrics csv file name to output
  N<sub>S</sub><sub>name</sub>              Susceptible array csv file name to output
  N<sub>P</sub><sub>name</sub>              Pre-symptomatic array csv file name to output
  N<sub>A</sub><sub>name</sub>              Asymptomatic array csv file name to output
  N<sub>Y</sub><sub>name</sub>              Symptomatic array csv file name to output
  N<sub>R</sub><sub>name</sub>              Recovered array csv file name to output
  time<sub>name</sub>             Time array csv file name to output

options:
  -h, &ndash;help            show this help message and exit
  -seed SEED            value used to initialize random number generator
  -N<sub>S0</sub> N<sub>S0</sub>            initial amount of susceptible population
  -N<sub>P0</sub> N<sub>P0</sub>            initial amount of pre-symptomatic population
  -N<sub>A0</sub> N<sub>A0</sub>            initial amount of asymptomatic population
  -N<sub>Y0</sub> N<sub>Y0</sub>            initial amount of symptomatic population
  -N<sub>R0</sub> N<sub>R0</sub>            initial amount of recoverd population
  -s S                  the fraction of contact that susceptible members will reduce
  -p P                  the fraction of contact that pre-symptomatic members will reduce
  -a A                  the fraction of contact that asymptomatic members will reduce
  -y Y                  the fraction of contact that symptomatic members will reduce
  -cycles CYCLES        the number of cycles the simulation will run
  -avg<sub>steps</sub> AVG<sub>STEPS</sub>  the number of equally distant in time averages we will be computing over the simulation time


<a id="org998d5b4"></a>

# Example Outputs

Ensure that you have followed the installation instructions beforehand:
Note: Outputs can be slightly different due to stochastic properties.
$ python simulation.py plot.png metrics.csv N<sub>S.csv</sub> N<sub>P.csv</sub> N<sub>A.csv</sub> N<sub>Y.csv</sub> N<sub>R.csv</sub> time<sub>name.csv</sub>
$ cat ./metrics.csv

<table border="2" cellspacing="0" cellpadding="6" rules="groups" frame="hsides">


<colgroup>
<col  class="org-left" />

<col  class="org-right" />

<col  class="org-right" />

<col  class="org-right" />

<col  class="org-right" />

<col  class="org-left" />
</colgroup>
<tbody>
<tr>
<td class="org-left">0</td>
<td class="org-right">1</td>
<td class="org-right">2</td>
<td class="org-right">3</td>
<td class="org-right">4</td>
<td class="org-left">5</td>
</tr>

<tr>
<td class="org-left">name</td>
<td class="org-right">mean</td>
<td class="org-right">std</td>
<td class="org-right">min</td>
<td class="org-right">max</td>
<td class="org-left">conf</td>
</tr>

<tr>
<td class="org-left">peak<sub>infections</sub></td>
<td class="org-right">1607.7</td>
<td class="org-right">54.329826062670215</td>
<td class="org-right">1512.0</td>
<td class="org-right">1721.0</td>
<td class="org-left">(np.float64(1592.1028752724171), np.float64(1623.297124727583))</td>
</tr>

<tr>
<td class="org-left">peak<sub>times</sub></td>
<td class="org-right">95548.3358400622</td>
<td class="org-right">9304.133210121096</td>
<td class="org-right">75886.2527655324</td>
<td class="org-right">121690.99321655945</td>
<td class="org-left">(np.float64(92877.28502548918), np.float64(98219.38665463521))</td>
</tr>

<tr>
<td class="org-left">attack<sub>rates</sub></td>
<td class="org-right">0.908452</td>
<td class="org-right">0.006205521412419749</td>
<td class="org-right">0.8948</td>
<td class="org-right">0.9214</td>
<td class="org-left">(np.float64(0.9066705054051609), np.float64(0.9102334945948392))</td>
</tr>
</tbody>
</table>

[Plot of default run](examples/plot.png)

$ python simulation.py skeptic<sub>plot.png</sub> skeptic<sub>metrics.csv</sub> skeptic<sub>N</sub><sub>S.csv</sub> skeptic<sub>N</sub><sub>P.csv</sub> skeptic<sub>N</sub><sub>A.csv</sub> skeptic<sub>N</sub><sub>Y.csv</sub> skeptic<sub>N</sub><sub>R.csv</sub> ignorance<sub>time</sub><sub>name.csv</sub> -s 0.80 -p 0.80 -a 0.80 -y 0.30
$ cat ./skeptic<sub>metrics.csv</sub>

<table border="2" cellspacing="0" cellpadding="6" rules="groups" frame="hsides">


<colgroup>
<col  class="org-left" />

<col  class="org-right" />

<col  class="org-right" />

<col  class="org-right" />

<col  class="org-right" />

<col  class="org-left" />
</colgroup>
<tbody>
<tr>
<td class="org-left">0</td>
<td class="org-right">1</td>
<td class="org-right">2</td>
<td class="org-right">3</td>
<td class="org-right">4</td>
<td class="org-left">5</td>
</tr>

<tr>
<td class="org-left">name</td>
<td class="org-right">mean</td>
<td class="org-right">std</td>
<td class="org-right">min</td>
<td class="org-right">max</td>
<td class="org-left">conf</td>
</tr>

<tr>
<td class="org-left">peak<sub>infections</sub></td>
<td class="org-right">618.1</td>
<td class="org-right">56.63717860204549</td>
<td class="org-right">487.0</td>
<td class="org-right">752.0</td>
<td class="org-left">(np.float64(601.8404754829234), np.float64(634.3595245170767))</td>
</tr>

<tr>
<td class="org-left">peak<sub>times</sub></td>
<td class="org-right">183043.46592964712</td>
<td class="org-right">26366.99888389911</td>
<td class="org-right">141785.27822713408</td>
<td class="org-right">254357.7722571604</td>
<td class="org-left">(np.float64(175473.9704961474), np.float64(190612.96136314684))</td>
</tr>

<tr>
<td class="org-left">attack<sub>rates</sub></td>
<td class="org-right">0.641204</td>
<td class="org-right">0.02227233225326885</td>
<td class="org-right">0.5922</td>
<td class="org-right">0.6768</td>
<td class="org-left">(np.float64(0.6348100103758165), np.float64(0.6475979896241835))</td>
</tr>
</tbody>
</table>

[Plot of skeptic run](examples/skeptic_plot.png)

$ python simulation.py infecpop<sub>plot.png</sub> infecpop<sub>metrics.csv</sub> infecpop<sub>N</sub><sub>S.csv</sub> infecpop<sub>N</sub><sub>P.csv</sub> infecpop<sub>N</sub><sub>A.csv</sub> infecpop<sub>N</sub><sub>Y.csv</sub> infecpop<sub>N</sub><sub>R.csv</sub> infecpop<sub>time</sub><sub>name.csv</sub> -N<sub>Y0</sub> 500
$ cat ./infecpop<sub>metrics.csv</sub>

<table border="2" cellspacing="0" cellpadding="6" rules="groups" frame="hsides">


<colgroup>
<col  class="org-left" />

<col  class="org-right" />

<col  class="org-right" />

<col  class="org-right" />

<col  class="org-right" />

<col  class="org-left" />
</colgroup>
<tbody>
<tr>
<td class="org-left">0</td>
<td class="org-right">1</td>
<td class="org-right">2</td>
<td class="org-right">3</td>
<td class="org-right">4</td>
<td class="org-left">5</td>
</tr>

<tr>
<td class="org-left">name</td>
<td class="org-right">mean</td>
<td class="org-right">std</td>
<td class="org-right">min</td>
<td class="org-right">max</td>
<td class="org-left">conf</td>
</tr>

<tr>
<td class="org-left">peak<sub>infections</sub></td>
<td class="org-right">1889.76</td>
<td class="org-right">52.01098345542026</td>
<td class="org-right">1788.0</td>
<td class="org-right">1993.0</td>
<td class="org-left">(np.float64(1874.828573655607), np.float64(1904.691426344393))</td>
</tr>

<tr>
<td class="org-left">peak<sub>times</sub></td>
<td class="org-right">39536.204943631106</td>
<td class="org-right">2214.2513855250227</td>
<td class="org-right">35217.0101841283</td>
<td class="org-right">44553.39182111935</td>
<td class="org-left">(np.float64(38900.53283604108), np.float64(40171.877051221134))</td>
</tr>

<tr>
<td class="org-left">attack<sub>rates</sub></td>
<td class="org-right">0.9083520000000002</td>
<td class="org-right">0.006593155238578862</td>
<td class="org-right">0.8948</td>
<td class="org-right">0.9204</td>
<td class="org-left">(np.float64(0.9064592226425721), np.float64(0.9102447773574283))</td>
</tr>
</tbody>
</table>

[Plot of infecpop run](examples/infecpop_plot.png)

