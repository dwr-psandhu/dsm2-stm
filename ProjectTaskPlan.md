# Introduction #

This document is the larger scope To Do List for the DSM2 Sediment and Transport Module (STM) project.

The document will be frequently modified by the project participants to
  * help plan the project direction
  * break project goals into components the subtasks of which can be tracked as [Issues](Issues.md). The issue tags that should be used to associate tasks with components are given below.


# Status #
Currently, the transport module is largely ready for single channels except for some details involving boundaries... even there it is fairly accurate and usable, it just has issues with formal convergence with small step sizes and complex active boundaries.

# Components and Details #
The main areas for STM Single Channel Model are shown below

## 1D Transport Solver ##
The development of a 1D transport solver is largely complete with the exception of some boundary management with operator splitting on problems with active boundaries.

Issues relating to the transport solver should be marked Component-Solver


## Sediment ##
Sediment is a separate module in STM. Within the transport code it takes the form of a source term, albeit a complex one. This job is currently the responsibility of the UCD group.

Issue tag: Component-Sediment

## Example Driver for Single Channel ##
There is interest in creating a standalone driver to go along with flume experiments. This could be based on DSM2 or very simple Fortran reads. The DSM2-compatible driver should be done from inside DSM2, not from this project.

Issue tag: Please do not tag things with this category until a preliminary design is ready.

## Tests ##
We have two categories of tests: unit tests for transport and sediment and a suite of convergence tests at the system level.

Test reports are generated automatically from code, but getting them into spreadsheets for presentations is a manual cut-and-paste operation. Hence these reports rapidly get out of date.

Issue tag: Component-Tests