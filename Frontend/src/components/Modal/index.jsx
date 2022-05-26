import { Dialog, DialogTitle } from '@mui/material';

/**
 * Modal
 * @param isOpen boolean state passed as prop, if false returns null
 * @param setOpen function to set boolean state of isOpen
 * @param children mui dialog title / dialog content / dialog actions children component
 */
export const Modal = ({
  isOpen, setOpen, children, sx, maxWidth='sm'
}) => (
  <Dialog
    open={isOpen}
    onClose={() => setOpen(false)}
    scroll="paper"
    aria-labelledby="scroll-dialog-title"
    aria-describedby="scroll-dialog-description"
    PaperProps={{
      sx: {
        borderRadius: '20px',
        p: '1em',
        ...sx,
      },
    }}
    maxWidth={maxWidth}
  >
    {children}
  </Dialog>
);

export const ModalTitle = ({ children }) => (
  <DialogTitle
    id="scroll-dialog-title"
    sx={{
      fontSize: '1.7em',
      fontWeight: 'bold',
      p: '0.5em',
    }}
  >
    {children}
  </DialogTitle>
);
