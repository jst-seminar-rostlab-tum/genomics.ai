import { Button, Typography } from "@mui/material"
import { Box } from "@mui/system"
import ReactDOM from "react-dom"
import { colors } from "shared/theme/colors"

/**
 * Modal
 * @param message the Modal bold title message 
 * @param isOpen boolean state passed as prop, if false returns null
 * @param onClose function for Close button onClick event
 * @param children react children component
 */
const Modal = ({ 
  message, isOpen, onClose, children
}) => {
  if(!isOpen) return null
  return ReactDOM.createPortal(
    <Box sx={{ overflow: "hidden" }}>
      <Box sx={{ 
        position: "absolute", top: "0", left: "0",
        width: "100vw", height: "100vw",
        backgroundColor: "rgba(0, 0, 0, 0.6)",
        zIndex: "20",
        overflow: "hidden"
      }}/>
      <Box
        sx={{
          position: "absolute",
          top: "50%",
          left: "50%",
          transform: "translate(-50%, -50%)",
          backgroundColor: "white",
          p: "2em",
          borderRadius: "20px",
          boxShadow: "0px 4px 6px rgba(0, 0, 0, 0.2), 0px 0px 1px rgba(0, 0, 0, 0.3)",
          zIndex: "100"
        }}
        >
        <Typography fontSize="2em" fontWeight="bold">{message}</Typography>
        <Box sx={{ height: "1px", width: "100%", backgroundColor: colors.neutral[400], margin: "0.2em 0 1em 0" }}/>
        {children}
        {/* TODO Replace with new buttons */}
        <Button onClick={onClose}>Close</Button>
      </Box>
    </Box>
  , document.body)
}

export default Modal