import { CloudDone, CloudUpload } from "@mui/icons-material"
import { Box, Typography } from "@mui/material"
import { useCallback } from "react"
import { useDropzone } from "react-dropzone"
import { colors } from "shared/theme/colors"

const BeforeDrag = () => {
  return (
    <Box sx={{ 
      display: "flex", 
      flexDirection: "column", 
      alignItems: "center", 
      gap: "1em", 
      color: colors.primary[800],
      position: "absolute",
      left: "50%",
      top: "50%",
      transform: "translate(-50%, -50%)",
      color: 'inherit'
    }}>
      <CloudUpload sx={{ fontSize: 80 }} />
      <Typography fontWeight="bold" fontSize="1.2em" sx={{ textAlign: "center" }}>Drag and drop or browse files</Typography>
    </Box>
  )
}

const DuringDrag = () => {
  return (
    <Box sx={{ 
      display: "flex", 
      flexDirection: "column", 
      alignItems: "center", 
      gap: "1em", 
      color: colors.primary[500], 
      position: "absolute",
      left: "50%",
      top: "50%",
      transform: "translate(-50%, -50%)"
    }}>
      <CloudDone style={{ fontSize: 80 }}/>
      <Typography fontWeight="bold" fontSize="1.2em" sx={{ textAlign: "center" }}>Drop files here</Typography>
    </Box>
  )
}

/**
 * FileUpload
 * @param width default value is 100% -> set it like "100px" or "8em"
 * @param height default value is 100% -> set it like "100px" or "8em"
 * @param handleFileUpload function that will handle all accepted files
 * @param rejectionHandler
 */
const FileUpload = ({ width = "100%", height = "100%", handleFileUpload, validator, rejectionHandler }) => {

  const onDropAccepted = useCallback( acceptedFiles => {
    handleFileUpload(acceptedFiles)
  }, [])

  const onDropRejected = useCallback(() => {
    rejectionHandler();
  }, [])

  const { getRootProps, getInputProps, isDragActive } = useDropzone({ onDropAccepted, onDropRejected, validator: validator })

  return (
    <Box
      {...getRootProps()}
      sx={{
        width: width,
        height: height,
        boxShadow: "0px 4px 6px rgba(0, 0, 0, 0.25), 0px 0px 1px rgba(0, 0, 0, 0.32)",
        p: "3em",
        borderRadius: "20px",
        cursor: "pointer",
        position: "relative",
        backgroundColor: colors.primary[100],
        '&:hover': {
          boxShadow: "inset 0px 8px 12px rgba(0, 0, 0, 0.25), 0px 0px 2px rgba(0, 0, 0, 0.32)",
          transition: '0.4s',
          backgroundColor: colors.primary[300],
          color: 'white'
        }
      }}
    >
     <input {...getInputProps()}/>
      {
        isDragActive ? 
        <DuringDrag />
        :
        <BeforeDrag />
      } 
    </Box>
  )
}

export default FileUpload