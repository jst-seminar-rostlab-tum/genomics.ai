import { Circle } from "@mui/icons-material"
import { Box, Button, Divider, Typography } from "@mui/material"
import { colors } from "shared/theme/colors"
import tum from 'assets/landing-illustrations/tum-logo.png';
import rostlab from 'assets/landing-illustrations/rostlab.png';
import helmholtz from 'assets/landing-illustrations/helmholtz.png';
import dnaImage from 'assets/dna.png';

const WindowiOS = () => {
  return (
    <Box
      sx={{
        position: "relative",
        width: "100%",
        height: "100%",
        margin: "auto",
        backgroundColor: "rgba(20, 20, 20, 0.6)",
        color: "white",
        p: "0.7em",
        borderRadius: "20px",
        boxShadow: "0px 0px 10px rgba(255,255,255, 0.15)",
        marginTop: "4em",
        marginBottom: "4em"
      }}
    >
      <Box
        sx={{
          position: "absolute",
          display: "flex",
          flexDirection: "row",
          gap: "10px",
          left: "3%",
        }}
      >
        <Circle sx={{ width: "15px", color: colors.error.main }}/>
        <Circle sx={{ width: "15px", color: colors.warning.main }}/>
        <Circle sx={{ width: "15px", color: colors.success.main }}/>
      </Box>
      <Typography textAlign="center">genomics.ai</Typography>
      <Box sx={{ p: "1em 4em 1em 4em", marginTop: "2em" }}>
        <Typography fontSize="3.3rem" fontWeight="bold">
          AI-Driven Cell Type Annotation
        </Typography>
        <Typography fontSize="2.2rem" fontWeight="semibold">
          No Code. Just Results.
        </Typography>
        <Typography fontSize="1.2rem" fontWeight="light" color="#8193B2">
          Genomics.ai helps you visualize all of your single-cell sequencing data in a fast and easy way using neural networks.
        </Typography>
        <Button disableRipple variant="contained" sx={{ margin: "2em 0 2em 0", backgroundColor: colors.primary[400], borderRadius: "20px" }}>Sign up for genomics.ai</Button>
        
        {/* DIVIDER */}
        <Box sx={{ marginTop: "3em", display: "flex", flexDirection: "row", justifyContent: "space-between", alignItems: "center" }}>
          <Box sx={{ height: "1px", width: "40%", backgroundColor: "rgba(255,255,255, 0.5)" }} />
          <Typography>Our partners</Typography>
          <Box sx={{ height: "1px", width: "40%", backgroundColor: "rgba(255,255,255, 0.5)" }} />
        </Box>
      </Box>
      {/* BUTTON PLACEHOLDER UNTIL WE HAVE A BUTTON DESIGN */}
      <Box
        sx={{
          display: "flex",
          flexDirection: "row",
          justifyContent: "space-between",
          alignItems: "center",
          width: "80%",
          margin: "auto",
          p: "1em"
        }}
      >
        <Box>
          <img src={tum} alt="TUM" style={{ width: "6em", height: "auto" }}/> 
        </Box>
        <Box>
          <img src={rostlab} alt="ROSTLAB" style={{ width: "10em", height: "auto" }}/> 
        </Box>
        <Box>
          <img src={helmholtz} alt="HELMHOLTZ" style={{ width: "6em", height: "auto" }}/> 
        </Box>
      </Box>
    </Box>
  )
}

export default WindowiOS