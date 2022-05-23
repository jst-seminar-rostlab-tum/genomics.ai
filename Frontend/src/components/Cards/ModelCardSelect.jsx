import { Box, Button, Container, Typography } from "@mui/material"
import ModelInfo from "components/GeneMapper/ModelInfo"
import { Modal } from "components/Modal"
import { useEffect, useRef, useState } from "react"
import { colors } from "shared/theme/colors"
import { LearnMoreModelComponent } from "views/Explore/LearnMoreModel"

// Outlined Button specific to Models and Atlases
/**
 * OutlinedButton
 * Specific to Model and Atlas Cards
 * @param content text content to be displayed
 * @param link button href link  
 */
export const OutlinedButtonSelect = ({ content, onSelect, disabled=false }) => {
 return (
  <Button
    variant="outlined"
    disableRipple
    onClick={onSelect}
    sx={{
      color: !disabled ? colors.primary[100] : "#FFFFFF" , borderWidth: "2px", 
      borderColor: !disabled ? colors.primary[100] : "#FFFFFF", 
      borderRadius: "1.2rem", justifyContent: 'center', textAlign: 'center',
      ":hover": { 
        borderColor: !disabled ? '#01579B' : "#696969", borderWidth: "2px", 
        backgroundColor: !disabled ? colors.primary[100] : "#FFFFFF", 
        color: !disabled ? '#01579B' : "#696969", transition: '0.5s' 
    }
    }}
  >
    {content}
   </Button>
 )
}

/**
 * ModelCard
 * @param width default value 100% of parent
 * @param height default value 100% of parent
 * @param title
 * @param description
 * @param mapLink onHover button Map url
 * @param learnMoreLink onHover button Learn More url
 */
export const ModelCardSelect = ({ 
  width = "100%", height = "100%", title, description, onSelect, 
  selected, learnMoreLink, modelObject={}, disabled=false
}) => {

  const [hover, setHover] = useState(false)
  const ref = useRef()
  const [flexDir, setFlexDir] = useState("column")
  const [modelInfoOpen, setModelInfoOpen] = useState(false);

  useEffect(() => {

    // checks if the parent element is wider than it is longer
    // if it is, converts the flex direction to row
    if(ref.current.clientWidth > ref.current.clientHeight) {
      setFlexDir("row")
    }
  }, [])

  return (
    <Box
      sx={{
        width: width,
        height: height
      }}
    >
      <Box
        ref={ref}
        sx={{ 
          width: "100%", 
          height: "100%", 
          position: "relative", 
          cursor: "pointer", 
          display: "flex", 
          flexDirection: "column", 
          justifyContent: "center"
        }}
        onMouseEnter={() => setHover(true)}
        onMouseLeave={() => setHover(false)}
        >  
        {
          hover && 
          <Box 
          style={{ background: !disabled ? "linear-gradient(#4F83CC, #01579B)" : "linear-gradient(#8A8A8A, #565656)" }} 
          sx={{ 
            position: "absolute", 
            width: "100%", 
            height: "100%", 
            borderRadius: "1.2rem",
            p:"1rem",
            opacity: 0.9
          }}>
            <Box
              sx={{
                display: "flex",
                flexDirection: flexDir,
                gap: "1rem",
                // styles below center the buttons to the very center of the component
                position: "absolute",
                top: "50%",
                left: "50%",
                transform: "translate(-50%, -50%)"
              }}  
              >
              {!disabled && <OutlinedButtonSelect content="Select" onSelect={() => onSelect(modelObject)}/>}
              <OutlinedButtonSelect content="Learn More" onSelect={() => setModelInfoOpen(true)} disabled={disabled}/>
            </Box>
          </Box>
        }
      <Box sx={{
          width: '100%',
          height: "100%",
          display: "flex",
          flexDirection: "column",
          p: "1.2rem",
          boxShadow: hover ? "none" : "0px 0px 10px rgba(0, 0, 0, 0.15)",
          borderRadius: "1.2rem",
          borderStyle:"solid",
          borderColor: selected ? "#008BF5" : 'transparent',
          borderWidth:"4px",
        }}
      >
        <Typography sx={{ fontSize: "1.4rem", fontWeight: "bold", color: !disabled ? "#000000" : "#8A8A8A"}}>{title}</Typography>
        <Typography 
          className="modelDescription" 
          sx={{ 
            fontSize: "1rem", 
            color: !disabled ? colors.neutral[800] : "#8A8A8A",
            display: '-webkit-box',
            overflow: 'hidden',
            WebkitBoxOrient: 'vertical',
            WebkitLineClamp: 2,
          }}>{description}</Typography>
      </Box>
      </Box>
      <ModelInfo id={modelObject._id} open={modelInfoOpen} setOpen={setModelInfoOpen} />
    </Box>
  )
}