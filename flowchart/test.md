mermaid
graph TB
software(VaspCZ) -->|Optimization and Static calculation|1.1
software(VaspCZ) -->|Transition state calculation| 2.4
software(VaspCZ) --> 3.1
subgraph Test module
3.1(ENCUT test) -->
3.2(KPOINTS test)
end
3.2-->|Results transfer|1.1
subgraph NEB module
2.1[Optimization-Static]
2.2{Static-NEB}
2.3[NEB-Vibration Analysis]
2.7(Check convergence)
2.4[Deal with input files]
2.5((Check NEB results))
2.6((Check NEB vib. results))l
2.1-->2.2
2.2-->2.3
2.4-->2.1
2.2-->|Get results|2.5
2.3-->|Get results|2.6
2.7-.->2.2
2.2-.->2.7
end
subgraph OS module
1.1[Deal with input files] -->
1.2[Pre-check and submit job]
1.3((Check results))
1.2-->|Get results|1.3
end
1.3-->|Results transfer|2.1


