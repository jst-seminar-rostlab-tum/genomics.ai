import { Box, Typography, Button} from "@mui/material";
import UploadIcon from '@mui/icons-material/Upload';
import AutorenewIcon from '@mui/icons-material/Autorenew';
import AssignmentIcon from '@mui/icons-material/Assignment';
import CheckCircleOutlineIcon from '@mui/icons-material/CheckCircleOutline';
import Navbar from "components/Navbar";
import WindowiOS from "components/WindowiOS";
import Footer from "components/Footer";
import { useCallback, useState, useEffect, useRef } from "react";
import { colors } from "shared/theme/colors";
import graphic1 from 'assets/landing-illustrations/science.png';
import graphic2 from 'assets/landing-illustrations/upload.png';
import graphic3 from 'assets/landing-illustrations/processing.png';
import graphic4 from 'assets/landing-illustrations/results.png';

const Home = ({ setUser }) => {
  const [isLoginFormVisible, setLoginFormVisible] = useState(false)
  const [isRegistrationFormVisible, setRegistrationFormVisible] = useState(false)

  // const howItWorksBoxRef = useRef()
  // const [processLineInfo, setProcessLineInfo] = useState({left: 0, top: 0, length: 0})

  // function updateProcessLineInfo(){
  //   console.log(howItWorksBoxRef)

  //   const iconSize = 16

  //   const updateBox = howItWorksBoxRef.current.children[1] 
  //   const startIcon = updateBox.children[0] 
  //   const startLeft = updateBox.offsetLeft + startIcon.offsetLeft + iconSize / 2
  //   const startTop = updateBox.offsetTop + startIcon.offsetTop + iconSize / 2

  //   const checkReusltBox = howItWorksBoxRef.current.children[3] 
  //   const endIcon = checkReusltBox.children[0] 
  //   const endTop = checkReusltBox.offsetTop + endIcon.offsetTop - iconSize / 2

  //   const newProcessLineInfo = {left: startLeft, top: startTop, length: endTop - startTop}

  //   console.log(newProcessLineInfo.left, newProcessLineInfo.top, newProcessLineInfo.length)

  //   if(processLineInfo.left!=newProcessLineInfo.left && processLineInfo.top!=newProcessLineInfo.top && processLineInfo.length!=newProcessLineInfo.length) setProcessLineInfo(newProcessLineInfo)
  // }

  // useEffect(()=>{
  //   // updateProcessLineInfo()

  //   window.addEventListener("resize", updateProcessLineInfo)
  //   return _=>window.removeEventListener("resize", updateProcessLineInfo)
  // })

  const onLoginClicked = useCallback(() => {
    setRegistrationFormVisible(false)
    setLoginFormVisible(true)
  }, [setLoginFormVisible])

  const onSignUpClicked = useCallback(() => {
    setLoginFormVisible(false);
    setRegistrationFormVisible(true);
  }, [setRegistrationFormVisible])
  
  const onRegistrationFormClosed = useCallback(() => {
    setRegistrationFormVisible(false);
  }, [setRegistrationFormVisible]);

  return (
    <Box style={{ overflow: "hidden" }} sx={{ position: "relative" }}>
      {/* STARTING PAGE */}
      <Box sx={{width: "100vw", bgcolor: colors.primary[800], position: "relative", paddingBottom: "4em"}}>
        {/* NAVBAR HERE */}
        <Navbar />
        {/* IOS WINDOW */}
        <WindowiOS />
      </Box>
      {/* the Eclipse */}
      <Box
        sx={{
          position: "relative",
          width: "140vw",
          height: "5vw",
          left: "-20vw", 
          top: "-2vw",
          backgroundColor: "white",
          borderRadius: "50%",
          zIndex: "0"
        }}
      />
      {/* CONTAINER */}
      <Box
        sx={{
          position: "relative",
          width: "100vw",
          paddingLeft: { xs: "10px", md: "15%" },
          paddingRight: { xs: "10px", md: "15%" },
          top: "-2.5vw"
        }}
      >
        
        {/* WHAT WE DO */}
        <Box sx={{ width: "100%" }}>
          <Typography width="100%" sx={{ textAlign: "center" }} fontSize="2em" fontWeight="bold">What we do</Typography>
          <Box sx={{ display: "flex", flexDirection: "row", justifyContent: "space-between", alignItems: "center", width: "80%", margin: "auto" }}>
            <img style={{ width: "25em" }} src={graphic1} alt="Science" />
            <Box sx={{ width: "40%", p: "1em"  }}>
              <Typography fontSize="1.2em" fontWeight="bold">genomics.ai</Typography>
              <Typography margin="2em 0 2em 0">We help you visualize all of your single-cell sequencing data in a fast and easy way with the help of neural networks.</Typography>
              {/* WILL BE REPLACED BY FUTURE BUTTONS */}
              <Button variant="contained">Explore</Button>
            </Box>
          </Box>
        </Box>
        {/* HOW IT WORKS */}
        <Box sx={{ 
          width: "100%", 
          marginTop: "5em", 
          backgroundColor: colors.primary[800],
          p: "1em",
          color: "white",
          borderRadius: "20px",
          boxShadow: "0px 0px 10px rgba(0, 0, 0, 0.1)",
          position: "relative"
          }}>
          <Typography sx={{ textAlign: "center" }} fontSize="2em" fontWeight="bold">How it works</Typography>

          {/* UPLOAD */}
          <Box sx={{ position: "relative", gap: "1em", display: "flex", flexDirection: {xs: "column", sm: "row", md: "row", lg: "row", xl: "row"}, justifyContent: "space-between", alignItems: "center", width: "90%", margin: "auto", padding: "10% 10% 5% 10%" }}>
            <Box sx={{position: "absolute", top: "calc(5% + 8px)", left: "7px", width: "2px", height: "calc(95% - 8px)", bgcolor: colors.neutral[700]}} />
            <Box sx={{width: "16px", height: "16px", position: "absolute", top: "5%", left: "0%", borderRadius: "10px", border: "2px solid white", bgcolor: colors.primary[800]}}></Box>
            <Box sx={{width: "30px", height: "30px", position: "absolute", top: "50%", left: "0%", transform: "translate(-7px, -15px)", bgcolor: colors.primary[800]}}>
              <UploadIcon sx={{width: "20px", height: "20px", position: "absolute", top: "50%", left: "50%", transform: "translate(-50%, -50%)"}}/>
            </Box>
            <Box sx={{ width: {xs: "100%", sm: "50%", md: "50%", lg: "50%", xl: "50%"} }}>
              <Typography fontSize="1.2em" fontWeight="bold">Upload</Typography>
              <Typography color={colors.neutral[500]}>Upload your files containing the single-cell sequencing data.</Typography>
            </Box>
            <Box sx={{ width: {xs: "100%", sm: "50%", md: "50%", lg: "50%", xl: "50%"}, backgroundColor: "white", borderRadius: "20px" }}>
              <img style={{ width: "100%" }} src={graphic2} alt="Upload" />
            </Box>
          </Box>  

          {/* PROCESSING */}
          <Box sx={{ position: "relative", gap: "1em", display: "flex", flexDirection: {xs: "column", sm: "row", md: "row", lg: "row", xl: "row"}, justifyContent: "space-between", alignItems: "center", width: "90%", margin: "auto", padding: "5% 10%" }}>
            <Box sx={{position: "absolute", top: "0%", left: "7px", width: "2px", height: "100%", bgcolor: colors.neutral[700]}} />
            <Box sx={{zIndex: "1", width: "30px", height: "30px", position: "absolute", top: "50%", left: "0%", transform: "translate(-7px, -15px)", bgcolor: colors.primary[800]}}>
              <AutorenewIcon sx={{width: "20px", height: "20px", position: "absolute", top: "50%", left: "50%", transform: "translate(-50%, -50%)"}}/>
            </Box>
            <Box sx={{ width: {xs: "100%", sm: "50%", md: "50%", lg: "50%", xl: "50%"} }}>
              <Typography fontSize="1.2em" fontWeight="bold">Processing</Typography>
              <Typography color={colors.neutral[500]}>Your input data is processed by our machine learning model which performs the cell-type classification according to your specification.</Typography>
            </Box>
            <Box sx={{ width: {xs: "100%", sm: "50%", md: "50%", lg: "50%", xl: "50%"}, backgroundColor: "white", borderRadius: "20px" }}>
              <img style={{ width: "100%" }} src={graphic3} alt="Processing" />
            </Box>
          </Box>

          {/* CHECK RESULT */}
          <Box sx={{ position: "relative", gap: "1em", display: "flex", flexDirection: {xs: "column", sm: "row", md: "row", lg: "row", xl: "row"}, justifyContent: "space-between", alignItems: "center", width: "90%", margin: "auto", padding: "5% 10% 10% 10%" }}>
            <Box sx={{position: "absolute", top: "0%", left: "7px", width: "2px", height: "calc(95% - 8px)", bgcolor: colors.neutral[700]}} />
            <Box sx={{zIndex: "1", width: "16px", height: "16px", position: "absolute", top: "95%", left: "0%", transform: "translate(-0%, -100%)", bgcolor: colors.primary[800]}}>
              <CheckCircleOutlineIcon sx={{width: "16px", height: "16px", color: colors.secondary1[400]}} />
            </Box>
            <Box sx={{zIndex: "1", width: "30px", height: "30px", position: "absolute", top: "50%", left: "0%", transform: "translate(-7px, -15px)", bgcolor: colors.primary[800]}}>
              <AssignmentIcon sx={{width: "20px", height: "20px", position: "absolute", top: "50%", left: "50%", transform: "translate(-50%, -50%)"}}/>
            </Box>
            <Box sx={{ width: {xs: "100%", sm: "50%", md: "50%", lg: "50%", xl: "50%"} }}>
              <Typography fontSize="1.2em" fontWeight="bold">Check Results</Typography>
              <Typography color={colors.neutral[500]}>Your input data is processed by our machine learning model which performs the cell-type classification according to your specification.After the algorithm has processed the data you can specify the project you want your cell-type data to be associated with and view the results.</Typography>
            </Box>
            <Box sx={{ width: {xs: "100%", sm: "50%", md: "50%", lg: "50%", xl: "50%"}, backgroundColor: "white", borderRadius: "20px" }}>
              <img style={{ width: "100%" }} src={graphic4} alt="Check Results" />
            </Box>
          </Box>
        </Box>
      </Box>

      {/* FOOTER */}
      <Footer />
    </Box>
  )
}

export default Home