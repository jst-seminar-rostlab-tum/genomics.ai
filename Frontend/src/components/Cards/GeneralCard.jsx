import { Box } from "@mui/material"

/**
 * GeneralCard
 * @param width default value 100%
 * @param height default value 
 * @param content text content
 * @param component if you provide a custom component it will show below content
 */
export const GeneralCard = ({ width = "100%", height = "100%", content, component } ) => {

  return (
    <Box
      sx={{
          width: width,
          height: height
      }}
    >
      <Box sx={{ boxShadow: "0px 0px 4px rgba(0,0,0, 0.25)", p: "0.8rem", borderRadius: "0.625rem", display: "flex", flexDirection: "row" }}>
        {content}
        {/* If the component is defined, show it to the screen */}
        {component && component}
      </Box>
    </Box>
  )
}