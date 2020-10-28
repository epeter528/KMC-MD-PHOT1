subroutine EQUIL_RELAXATION

           call system('rm -f equil.tpr')
           call system('/loctmp/pee18323/GRO_bin/bin/grompp -f fullmd_solequil.mdp -c md_cheater2.gro -p LOV2gen.top -o equil.tpr')
           call system('mpirun -np 4 /loctmp/pee18323/GRO_bin/bin/mdrun -v -s equil.tpr -c start.gro -nice 0')
           call system('rm -f #*#')

end 