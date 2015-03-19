# Introduction #

This module of DSM2 is a transport code modeling advection, diffusion, sediment processes and reactions.

# Recommended Work Plan #

1. Get doxygen working. The documentation is automatically built by the MS Visual Studio solution as long as doxygen is on path. See below on Visual Studio.<br><br>
2. Code a function-test pair and add it to test_transport_driver(). I suggest you either<br>
(1) adjust your diffusion code stylistically and write a matching unit test (see below on unit tests) or<br>
(2) complete the extrapolate() function in the advect module and write a test for that. See below for how to write good unit tests.<br><br>
3. Move the advection.replace_boundary_flux() subroutine out of the advection model since it it is driver/application specific. This step will teach you about interfaces and callbacks. Discuss the design with Nicky and Kevin and Eli.<br><br>
4. After that, start looking at which functions in the advection module are not really implemented. Work in function-test pairs (they often evolve together). As much as you can, put "asserts" in your testing function that are needed regardless of whether you can pass them. The Tests, and in some ways the documentation, are our way of communicating expected behavior.<br><br>
5. You will be at the same point you are in Matlab when you code the uniform flow advection test. However, we will not want you to code this or turn it in until the other testing (extrapolate(), conservative update) are well tested.<br><br>
6. Do not significantly change the granularity or form of these functions without a design meeting. It is recommended that you discuss the API (names and arguments) of routines and the basic flow as a skeleton before you code along with the tests. If you pass the tests and have a functional form we can use, you have succeeded at the task by definition. In the meantime:<br><br>
7. I recommend you work in this order, but you may rearrange the diffusion part as you wish:<br>
(1) code and test existing skeleton routines in the transport library<br>
(2) code test for uniform flow advection no source or diffusion<br>
-- achieve comfort with convergence testing<br>
(3) code test for non-uniform flow advection (handle flow/area consistency issues)<br>
(4) add linear source and test convergence<br>
(5) code test for diffusion<br>
-- explicit operator<br>
-- implicit solver<br>
------ convergence test<br>
(6) code test for advective-diffusion<br>
(7) meet on channel network challenges<br>
(8) modify for sediment<br>
(9) complete channel network stuff<br>

<h1>Visual Studio</h1>

Hopefully we will just give you a complete project that has most of the settings right. To get this to compile you need to make sure you have:<br>

1. Doxygen installed and on path. To get the doc project to generate, you need to have doxygen installed. Make sure it is on path (test by typing doxygen at the command line). You can add to the path on your computer or set it up in Visual Studio's Tools > Options > Projects and Solutions > VC++ Di drectories (for executables) Note the doc project has a custom build step (cd's to /doc and does documentation.bat) output for the custom build step is /doc/html/index.html You can view /doc/html/index.html in Visual Studio by right clicking it and choosing "View in Browser"<br><br>
2. Set up the dependencies (hopefully we will do this for you)<br>
---- test_transport on transport and fruit<br>
---- sandbox_application on transport (you can alter as desired)<br><br>
3. Set the fortran "additional include directories" in sandbox_application and test_transport_driver to pick up the fruit and transport module files. Hopefully we will do this for you.<br><br>
4. The sandbox is like scratch paper – it is an executable project to do whatever you want.<br><br>

<h1>Module Divisions: Library versus Application</h1>

You will notice that the main computational code is in a static library called /transport. In the future we may have other static library modules called /sediment /nonconservative and to deal with I/O.<br><br>

By and large, the drivers (test code and application) know about the libraries and not vice versa. The drivers do the setup and request the transport library to solve individual time steps. There are a couple nuanced points about this:<br><br>

1. The state_variables module properly allocates and houses state variables. It is like a utility for drivers tests and applications. You should use the state_variables module and then pass the variables into the solver as arguments.<br><br>
2. The advection.advect() function takes a callback for a source term that must be provided by the driver/application. The reason I did this is that we want to be able to stipulate different source terms for different tests. We will probably want to implement advection.replace_boundary_flux() the same way.<br><br>
3. When you use callbacks, you should do so at the level of the whole grid. In other words compute_source is perfectly efficient as long as its job is to compute source over the grid, not a point.<br><br>
4. Other common jobs done by the driver do not affect the interior of advect(). For instance, the initial condition, and hydro_data(). These have to be user supplied. I have tried to include interface blocks for these routines in the library, but not an implementation.<br><br>

<h1>Why All The <i>lo and</i>hi Data Structures??? Isn't this Wasteful?</h1>

You may notice that I have included a lot of data structures that look like flux_lo(icell) and flux_hi(icell) and in some cases (but not all) you may feel that there is a redundancy because it seems like flux_hi(icell) == flux_hi(icell+1). This, in fact will be true for some data structures for some data. However, when we get to channel networks it will not be true. At that point we will be stuck with a choice between:<br><br>

-- face-centered arrays whose mappings are non-obvious, ie: flux(hi_face(ncell)).<br>
-- flux_hi(ncell), which is still perfectly well-defined still.<br><br>

<h1>Tests</h1>

I put the Fruit testing framework here directly as a static library to reduce complexity of the project files. We can move it later.<br><br>

The test_transport project is the first battery of tests. The tests in test_gradient.f90 are fairly complete for that module. You will find a test driver called test_transport_driver.f90 The driver program there is where you would add new tests. Try to have one for every medium sized routine. For algorithms, the best integration tests are either accuracy compared to a known solution or (better) convergence.<br><br>

Remember, unit tests are silent. They should also be complete and catch "corner" cases. Look at the flux limiter test. I spot checked one "ordinary" case in the middle of the array in a gentle area where the limited flux is just the centered difference. More importantly, I tested both ends of the array and I looked at reversals of sign in both directions and places where the limiter is in use due to big jumps on the hi and lo side.<br><br>

This is usually a sobering experience.<br><br>

<h1>Todo Tags</h1>

In the code you should use !todo: for any notes to yourself about things to do later. You will already see a few in there. Do not use any other format for this job.<br><br>

<h1>Style Notes:</h1>

Please try to adhere to these guidelines as well as trying to get a sense of the code already provided and matching its style. Although we encourage you to "not worry about it" when you turn in code every other week, it would be nice if you took 1-2 hours to go over this as a checklist.<br><br>

-- Try to follow the subroutine style. The gradient module is probably the best example for routines and state_variables module shows how to encapsulate variables and document them.<br><br>
-- All routines are implicit none. Write this out, do not use compiler options.<br><br>
-- Do not pass data to subroutines by "use" statements. Put data in the signature. Imagine that you want to test the routine later out of context. Long signatures are fine.<br><br>
-- Routine names are lower_underscore. Try to name them with one well chosen word or two, avoiding redundancy. Do not abbreviate unless:<br>
------ the abbreviation is very easy to interpret and the original is taxing. For instance, conc is fine for concentration.<br>
------ the abbreviation is standard English, in which case you should always use it: e.g. abbr for abbreviation.<br>
------ the abbreviation is well established in the code. For instance I use lo and hi to show that a quantity is on face on the lo (n-1) or hi (n+1) side of the cell.<br><br>
-- Variables are also lower_underscore.<br><br>
-- Some constants associated with precision like LARGEREAL and GTM_REAL are all caps.<br>
------ OK, this is a bit consistent but the typesafe scalar values for numbers (one, two, half, etc) in the gtm_precision module are lower because... well it just looks way better<br><br>
-- Use the typesafe scalar values for numbers that have them: x = two*y. This avoids having to keep track of the "d" for double precision and prevents real-double conversion.<br><br>
-- Prefer subroutines to functions.<br><br>
-- Argument lists should not be long. Use continutation lines (&) when there is more than one line full. See the examples.<br><br>
-- Section off the arguments to all routines. See gradient.f90 for examples<br><br>
-- Section off routines from one another using //////. See gradient.f90<br><br>
-- Declare the intent for every variable.<br><br>
-- Use real(GTM_REAL) for floating point. Avoid hardwiring.<br><br>
-- Do not add things like "last modified" that have to be maintained by hand. Let version control do this.<br><br>
-- Put the licence at the top.<br><br>
-- Initialize variables to LARGEREAL or put LARGEREAL in indexes that won't get touched. For instance, if the derivative and value arrays are the same side and you do a "lo-side" difference, there is nothing to put in the first index. Set that to LARGEREAL. The reason for this is that if it gets initialized to something reasonable (the compiler often choses zero, at least in debug mode) it can lead to bugs that look deceptively reasonable.<br><br>
-- In some cases you may see a_prepended to an argument name. It means "argument". This was done when a global or module variable clashes with the argument name. It is not required (local names take precedence) but it prevents ambiguity.<br><br>
-- Indent with spaces, and consistently within a file. Use F90 free-form syntax. Again, see gradient.f90.<br><br>
-- Try not to check in code with commented out sections.<br><br>
Avoid comments when clear code names will suffice. Comments should outline the intent of little blocks of code when that intent is not obvious.