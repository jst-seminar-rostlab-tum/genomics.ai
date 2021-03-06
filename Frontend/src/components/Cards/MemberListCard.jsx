/* eslint-disable */

import { Box } from "@mui/material"

/**
 * GeneralCard
 * @param width default value 100%
 * @param height default value 
 * @param flex default value "column"
 * @param bg background color default value "white"
 * @param hoverStyles similar to sx prop, passed to ":hover"
 * 
 * Usage: 
 *  <GeneralCard>
 *    Your components
 *  </GeneralCard>
 */
export const MemberListCard = ({ width = "100%", height = "100%", children, flex = "column", bg = "white", hoverStyles = {}, padding = "0.5em", marginLeft = "0.2em" }) => {

  return (
    <Box
      sx={{
        width: width,
        height: height,
        backgroundColor: bg,
        borderRadius: "0.625rem",
        ":hover": hoverStyles
      }}
    >
      <Box sx={{ boxShadow: "0px 0px 4px rgba(0,0,0, 0.25)", p: padding, borderRadius: "0.625rem", display: "flex", flexDirection: flex }}>
        {children}
      </Box>
    </Box>
  )
}