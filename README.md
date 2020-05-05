# MarkovChain
The implementation of Markov Chain and HMM.

# Prepareation
- Create a folder named 'data', then put the file into this folder.
- In main file, you need to specify the following arguments.
    -data           : src path
    --entryname     : entry name
    -startpattern   : Starting with
    -endpattern     : Ending with
    -from           : starting from 'from' line index
    -to             : ending to 'to' line index
- NOTE # if the startpattern and endpattern is not found, this program will segmentation fault.


example:
```
> ./main --data ./data/NC_000006.12_chromosome_6.txt \
--entryname ">NC_000006.12 Homo sapiens chromosome 6, GRCh38.p13 Primary Assembly" \
--startpattern ttggtaccat \
--endpattern CTTTGCCTG \
--from 100000 \
--to 199999
```


# Illustration

- Note: Because I didn't find sequence described in the statements, I find another sequence which length is also approximately 100k long, starting with 'ttggtaccat' and ending with 'CTTTGCCTG',


## Simple Markov Chain

- Zero order markov chain
    - Log 2 Probability : -295341.059760


- A, T, C, G sta probability\

| a    | t    | c    | g   |
| ---- | ---- | ---- | --- |
| 2.929685e-01 | 3.192353e-01 | 1.970112e-01 |  1.907850e-01   |


- First order markov chain
    - Log 2 Probability : -584165.254544

```
                   a            t            c            g
    a   9.635000e-02 8.400391e-02 4.858698e-02 6.402956e-02 
    t   7.135216e-02 1.107892e-01 6.298632e-02 7.410976e-02 
    c   7.042188e-02 7.183058e-02 4.756367e-02 7.189703e-03 
    g   5.484640e-02 5.260710e-02 3.787552e-02 4.545726e-02 
```






- Second order markov chain
    - Log 2 Probability : -1158744.584939
```
                      aa          at          ac          ag          ta          tt          tc          tg          ca          ct          cc          cg          ga          gt          gc          gg
          aa  1.4120e-02  1.0559e-02  5.5085e-03  6.7710e-03  7.7744e-03  8.4655e-03  4.5185e-03  6.4255e-03  5.6746e-03  4.8839e-03  2.8373e-03  5.9803e-04  6.0002e-03  5.0235e-03  3.1563e-03  4.0533e-03
          at  6.9970e-03  7.8741e-03  3.5018e-03  3.7277e-03  6.6514e-03  1.1841e-02  5.1763e-03  4.3058e-03  5.1231e-03  5.3491e-03  3.0898e-03  3.8540e-04  5.8142e-03  5.9072e-03  3.7742e-03  4.5118e-03
          ac  5.2427e-03  5.7677e-03  4.0134e-03  4.7045e-03  3.3888e-03  5.2826e-03  3.2692e-03  4.4454e-03  3.9204e-03  3.5882e-03  2.4453e-03  3.1895e-04  5.4487e-04  7.4422e-04  4.5185e-04  4.5849e-04
          ag  7.1631e-03  4.7510e-03  3.2094e-03  4.6646e-03  3.7277e-03  5.3823e-03  3.1363e-03  4.4188e-03  4.4121e-03  4.1596e-03  3.1098e-03  3.9869e-04  4.9171e-03  3.8673e-03  3.6347e-03  3.0699e-03
          ta  8.9439e-03  6.9970e-03  2.9636e-03  3.4952e-03  6.5185e-03  8.5053e-03  3.4686e-03  4.7842e-03  4.7776e-03  3.9736e-03  2.4785e-03  4.7178e-04  4.3058e-03  3.8606e-03  2.5582e-03  3.2493e-03
          tt  7.7810e-03  7.7279e-03  3.8407e-03  4.6248e-03  9.6748e-03  1.6333e-02  8.9572e-03  8.1332e-03  6.9239e-03  9.5486e-03  6.2461e-03  4.0533e-04  5.4354e-03  6.4255e-03  4.0865e-03  4.6314e-03
          tc  5.3424e-03  5.9471e-03  4.2194e-03  4.9769e-03  4.5583e-03  8.0734e-03  5.8408e-03  5.9604e-03  5.5550e-03  6.2660e-03  4.1663e-03  4.2527e-04  4.9836e-04  5.5816e-04  3.0566e-04  3.3888e-04
          tg  6.3723e-03  5.2029e-03  3.1363e-03  5.0500e-03  5.1497e-03  7.3358e-03  4.0467e-03  5.3291e-03  4.9437e-03  5.3092e-03  4.1131e-03  4.2527e-04  4.7377e-03  4.7842e-03  3.4752e-03  4.6580e-03
          ca  6.8973e-03  4.8440e-03  2.9237e-03  3.7942e-03  4.6115e-03  6.7445e-03  3.6879e-03  5.1098e-03  5.7411e-03  4.6248e-03  3.1430e-03  7.2428e-04  4.8374e-03  4.5118e-03  3.9869e-03  4.2925e-03
          ct  4.0932e-03  3.8806e-03  2.4453e-03  2.8971e-03  4.2859e-03  8.3060e-03  5.6082e-03  4.6513e-03  5.1630e-03  5.5816e-03  4.7377e-03  6.1132e-04  4.7710e-03  5.5351e-03  4.2660e-03  4.9969e-03
          cc  4.3058e-03  4.9570e-03  3.5815e-03  4.5716e-03  2.8440e-03  5.3291e-03  4.3656e-03  4.8640e-03  4.3856e-03  3.6746e-03  2.3456e-03  5.5152e-04  4.2527e-04  4.5849e-04  3.9204e-04  5.0500e-04
          cg  4.3856e-04  4.2527e-04  3.3224e-04  5.7145e-04  4.7842e-04  6.1796e-04  4.3856e-04  7.1099e-04  3.5882e-04  4.8507e-04  5.1165e-04  1.3290e-04  3.8540e-04  4.5185e-04  3.8540e-04  4.6513e-04
          ga  6.9970e-03  4.7842e-03  2.5981e-03  4.1729e-03  3.1961e-03  4.2593e-03  2.2725e-03  3.6879e-03  3.5350e-03  2.9038e-03  1.8140e-03  4.0533e-04  4.6447e-03  3.2692e-03  2.3788e-03  3.8938e-03
          gt  3.5284e-03  3.7942e-03  1.9137e-03  2.7244e-03  3.3623e-03  6.6182e-03  3.3822e-03  3.4819e-03  3.2759e-03  3.9603e-03  2.3390e-03  2.9902e-04  3.7410e-03  3.9935e-03  2.6646e-03  3.5151e-03
          gc  3.5682e-03  3.4819e-03  2.4187e-03  3.3755e-03  2.5250e-03  4.1663e-03  2.6180e-03  4.2992e-03  3.5550e-03  3.8739e-03  2.0001e-03  4.8507e-04  2.9902e-04  4.8507e-04  3.3888e-04  3.8540e-04
          gg  4.5783e-03  3.0367e-03  1.9801e-03  3.9005e-03  2.6048e-03  3.5084e-03  2.2526e-03  3.4553e-03  3.1297e-03  3.6546e-03  2.1795e-03  5.5152e-04  3.4553e-03  2.7177e-03  2.0200e-03  2.4054e-03
```
## Hidden Markov Chain

### Sequence information

- Starting with "ttggtaccat"
- Ending with "CTTTGCCTG"
- Length : ~ 100k base

### HMM arch

![Image of HMM](https://i.imgur.com/jTk7Uzk.png)

![Image of HMM](https://i.imgur.com/0yY3v0o.png)

### Log probability

- Forward Algorithm
    - Log 2 Probability : -3.159066e+05

- Backward Algorithm
    - Log 2 Probability : -3.159073e+05


- Applying learning algorithm
    - Log 2 Probability : -3.009859e+05

- State sequence
The most likely state sequence is written in the state_sequence.txt


### Another 100kb section
- Information
    - Start : "gagatccttc"
    - End   : "ctttaaaagaaaa"
- Forward Algorithm
    - Log 2 probability : -2.324345e+06

