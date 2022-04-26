import { Box, Typography } from "@mui/material"
import { useCallback } from "react"
import { useDropzone } from "react-dropzone"
import { CloudUpload, CloudDone } from "@mui/icons-material"
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
      transform: "translate(-50%, -50%)"
    }}>
      <CloudUpload style={{ fontSize: 80 }}/>
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
 */
const FileUpload = ({ width = "100%", height = "100%", handleFileUpload }) => {

  const onDrop = useCallback( acceptedFiles => {
    handleFileUpload(acceptedFiles)
  }, [])

  const { getRootProps, getInputProps, isDragActive } = useDropzone({ onDrop })

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
        position: "relative"
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