import { Box, Typography, Divider } from "@mui/material"
import { colors } from "shared/theme/colors"

/**
 * Usage:  
 *  <GeneralCardHeader>{your text}</GeneralCardHeader>
 */
export const GeneralCardHeader = ({ children }) => {
  return (
    <Box sx={{ width: "100%" }}>
      <Typography color={colors.neutral[800]} fontSize="1.5em" fontWeight="bold">{children}</Typography>
      <Divider orientation="horizontal" sx={{ borderBottomWidth: 2 }} />
    </Box>
  )
}