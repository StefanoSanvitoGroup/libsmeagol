      module parallel
C   A basic version of the parallel module with stripped mesh routines.
C   Routines to convert between local and global orbital indices were
C   simplified by Sergey Chulkov (2024) to minimize the number of
C   floating-point operations.
C
C   This module is based on SIESTA's parallel.f module written by
C   Julian Gale in October 1998.

      implicit none
      private

C   ScaLAPACK blocking factor.
C   Setting this value is a compromise between dividing up the orbitals
C   over the processors evenly to achieve load balancing and making the
C   local work efficient. Typically a value of about 10 is good, but
C   optimisation may be worthwhile. A value of 1 is very bad for any
C   number of processors and a large value may also be less than ideal.
      integer, save :: BlockSize = 8

      public :: SetBlockSize
      public :: GetNodeOrbs
      public :: GlobalToLocalOrb
      public :: LocalToGlobalOrb
      public :: WhichNodeOrb

      contains

      subroutine SetBlockSize(n)
C
C  Set ScaLAPACK blocking factor.
C  This routine is not part of the original SIESTA module.
C
      integer, intent(in) :: n

      BlockSize = n
      end subroutine SetBlockSize

      pure subroutine GetNodeOrbs(NOrb, Node, Nodes, NOrbNode)
C
C  Calculates the number of orbitals stored on the local Node.
C
C  Input :
C
C  integer NOrb     = The total number of orbitals in the calculation
C  integer Node     = The local processor
C  integer Nodes    = The total number of processors
C
C  Output :
C
C  integer NOrbNode = The number of orbitals stored on this Node - if zero
C                     on input then calculated otherwise left unchanged

C Passed arguments
      integer, intent(in)  :: NOrb, Node, Nodes
      integer, intent(out) :: NOrbNode

C Local variables
      integer :: FullBlocks, NBlocks, RemainedBlocks, RemainedOrb

C-----------------------------
C  Block-cyclic distribution -
C-----------------------------

C  Total number of blocks across all nodes
      NBlocks = NOrb/BlockSize

C  Number of orbitals in the last partially-filled block
      RemainedOrb = NOrb - NBlocks*BlockSize

C  Number of full blocks on all nodes
      FullBlocks = NBlocks/Nodes

C  Remaining number of complete clocks
      RemainedBlocks = NBlocks - FullBlocks*Nodes

C  Local number of orbitals
      NOrbNode = FullBlocks*Nodes*BlockSize +
     .           MERGE(BlockSize, 0, Node < RemainedBlocks) +
     .           MERGE(RemainedOrb, 0, Node == RemainedBlocks)

      end subroutine GetNodeOrbs

      pure subroutine GlobalToLocalOrb(GOrb, Node, Nodes, LOrb)
C
C  Converts the global index of an orbital to its local index on a Node.
C  Returns 0 if the orbital is not local to the Node.
C
C  Input :
C
C  integer GOrb   = global orbital index
C  integer Node   = local processor number
C  integer Nodes  = global number of processors
C
C  Output :
C
C  integer LOrb   = local orbital index

C Passed arguments
      integer, intent(in)  :: GOrb, Node, Nodes
      integer, intent(out) :: LOrb

C  Local variables
      integer :: GBlock, FullBlocks

C-----------------------------
C  Block-cyclic distribution -
C-----------------------------

C  Global block number
      GBlock = (GOrb-1)/BlockSize

C  Number of full blocks on all nodes
      FullBlocks = GBlock/Nodes

      if (GBlock-FullBlocks*Nodes == Node) then
C        MOD(GBlock, Nodes) == Node : the local orbital

         LOrb = GOrb-BlockSize*(GBlock-FullBlocks)
      else
C        The orbital is local to some other node

         LOrb = 0
      end if
      end subroutine GlobalToLocalOrb

      pure subroutine LocalToGlobalOrb(LOrb, Node, Nodes, GOrb)
C
C  Converts the local index of an orbital on a Node to its global index.
C
C  Input :
C
C  integer LOrb   = local orbital index
C  integer Node   = local processor number
C  integer Nodes  = global number of processors
C
C  Output :
C
C  integer GOrb   = global orbital index

C Passed arguments
      integer, intent(in)  :: LOrb, Node, Nodes
      integer, intent(out) :: GOrb

C  Local variables
      integer :: LBlock

C  Local block number
      LBlock = (LOrb-1)/BlockSize

C  Global orbital index
      GOrb = (LBlock*(Nodes-1) + Node)*BlockSize + LOrb

      end subroutine LocalToGlobalOrb

      pure subroutine WhichNodeOrb(GOrb, Nodes, Node)
C
C  Returns the Node number local to an orbital with
C  the global index GOrb.
C
C  Input :
C
C  integer GOrb   = global orbital index
C  integer Nodes  = total number of Nodes
C
C  Output :
C
C  integer Node   = Node where this orbital is stored locally

C Passed arguments
      integer, intent(in)  :: GOrb, Nodes
      integer, intent(out) :: Node

C  Local variables
      integer :: GBlock

C-----------------------------
C  Block-cyclic distribution -
C-----------------------------
C  Find global block number
      GBlock = (GOrb-1)/BlockSize

C  Find the Node number that has this block
      Node = MOD(GBlock,Nodes)
      end subroutine WhichNodeOrb

      end module parallel
