import { Dialog, DialogTitle } from "@mui/material"

/**
 * Modal
 * @param isOpen boolean state passed as prop, if false returns null
 * @param children mui dialog title / dialog content / dialog actions children component
 */
const Modal = ({
  isOpen, children
}) => {
  return (
    <Dialog
      open={isOpen}
      onClose={() => setOpen(false)}
      scroll="paper"
      aria-labelledby="scroll-dialog-title"
      aria-describedby="scroll-dialog-description"
      PaperProps={{
        sx: {
          borderRadius: "20px",
          p: "1em"
        }
      }}
    >
      {children} 
    </Dialog>
  )
}

export const ModalTitle = ({ children }) => {
  return (
    <DialogTitle id="scroll-dialog-title"
      sx={{
        fontSize: "1.7em",
        fontWeight: "bold",
        p: "0.5em"
      }}
    >
      {children}
    </DialogTitle>
  )
}

export default Modal