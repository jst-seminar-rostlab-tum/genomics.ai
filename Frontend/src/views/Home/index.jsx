/* eslint-disable */

import { Box, Typography, Button, Stack, TextField, TextareaAutosize } from "@mui/material";
import UploadIcon from '@mui/icons-material/Upload';
import AutorenewIcon from '@mui/icons-material/Autorenew';
import AssignmentIcon from '@mui/icons-material/Assignment';
import CheckCircleOutlineIcon from '@mui/icons-material/CheckCircleOutline';
import Navbar from "components/NavBar";
import WindowiOS from "components/WindowiOS";
import Footer from "components/Footer";
import LoginForm from 'components/LoginForm'
import RegistrationForm from 'components/RegistrationForm'
import { useState, useEffect, useRef, useContext } from "react";
import { colors } from "shared/theme/colors";
import graphic1 from 'assets/landing-illustrations/science.png';
import graphic2 from 'assets/landing-illustrations/upload.png';
import graphic3 from 'assets/landing-illustrations/processing.png';
import graphic4 from 'assets/landing-illustrations/results.png';
import CustomButton from "components/CustomButton";
import Input from 'components/Input/Input'
import { useHistory, useLocation } from "react-router-dom";
import PasswordForgetForm from "components/PasswordForgetForm";
import ContactForm from 'components/ContactForm';
import { LoginContext } from "shared/context/loginContext";

const Home = () => {

  const history = useHistory();
  const location = useLocation();

  const context = useContext(LoginContext)

  const onLoginClicked = () => {
    context.switchRegister(false)
    context.switchLogin(true)
  }

  const onSignUpClicked = () => {
    context.switchLogin(false);
    context.switchRegister(true);
  }

  //here we get the ref of the contact us, in order to be able to scroll to it
  const contactUsBoxRef = useRef()
  const executeScroll = () => contactUsBoxRef.current.scrollIntoView()

  useEffect(() => {
    if (location.state && location.state.contact_us) contactUsBoxRef.current.scrollIntoView()
  }, [])

  //we store the actual height of the Navbar, since we set the Navbar's position to fixed
  //it jumps out from the document flow => the height collapse
  //we need to reset it
  const [navbarHeight, setNavbarHeight] = useState(0)

  return (
    <Box style={{ overflow: "hidden" }} sx={{ position: "relative" }}>
      {context.loginVisible && <LoginForm />}
      {context.registerVisible && <RegistrationForm  />}
      {context.forgetVisible && <PasswordForgetForm />}
      {/* STARTING PAGE */}
      <Box sx={{ width: window.width, bgcolor: colors.primary[800], position: "relative", paddingBottom: "4em" }}>
        {/* NAVBAR HERE */}
        {/* the Box that contains the Navbar will collapse, so we reset the height */}
        <Box sx={{ height: navbarHeight, position: "relative" }}>
          <Navbar position="fixed" setNavbarHeight={setNavbarHeight} onLoginClicked={onLoginClicked} onSignUpClicked={onSignUpClicked} executeScroll={executeScroll} />
        </Box>
        {/* IOS WINDOW */}
        <WindowiOS onSignUpClicked={onSignUpClicked} />
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
          padding: { xs: "0% 5%", sm: "0% 5%", md: "0% 15%", lg: "0% 15%", xl: "0% 15%" },
          top: "-2.5vw"
        }}
      >

        {/* WHAT WE DO */}
        <Box sx={{ width: "100%", position: "relative", margin: "auto" }}>
          <Typography sx={{ textAlign: "center" }} fontSize="2em" fontWeight="bold">What we do</Typography>
          <Box sx={{ display: "flex", flexDirection: { xs: "column", sm: "row", md: "row", lg: "row", xl: "row" }, justifyContent: "space-between", alignItems: "center", width: "80%", gap: "1em", margin: "auto" }}>
            <Box sx={{ width: { xs: "100%", sm: "50%", md: "50%", lg: "50%", xl: "50%" }, backgroundColor: "white", borderRadius: "20px" }}>
              <img style={{ width: "100%" }} src={graphic1} alt="Science" />
            </Box>
            <Box sx={{ width: { xs: "100%", sm: "50%", md: "50%", lg: "50%", xl: "50%" } }}>
              <Typography fontSize="1.2em" fontWeight="bold">Query to reference dataset mapping using Transfer Learning</Typography>
              <Typography margin="2em 0 2em 0">ArchMap enables mapping by using a novel deep learning strategy called single-cell architectural surgery (ScArches).</Typography>
              <CustomButton onClick={() => history.push('/explore')}>Explore</CustomButton>
            </Box>
          </Box>
        </Box>
        {/* HOW IT WORKS */}
        <Box sx={{
          width: "100%",
          marginTop: "5em",
          backgroundColor: colors.primary[800],
          p: "1em 1em 5em 1em",
          color: "white",
          borderRadius: "20px",
          boxShadow: "0px 0px 10px rgba(0, 0, 0, 0.1)",
          position: "relative"
        }}>
          <Typography sx={{ textAlign: "center" }} fontSize="2em" fontWeight="bold">How it works</Typography>

          {/* UPLOAD */}
          <Box sx={{ position: "relative", gap: "1em", display: "flex", flexDirection: { xs: "column", sm: "row", md: "row", lg: "row", xl: "row" }, justifyContent: "space-between", alignItems: "center", width: "90%", margin: "auto", padding: "10% 10% 5% 10%" }}>
            <Box sx={{ position: "absolute", top: "calc(5% + 8px)", left: "7px", width: "2px", height: "calc(95% - 8px)", bgcolor: colors.neutral[700] }} />
            <Box sx={{ width: "16px", height: "16px", position: "absolute", top: "5%", left: "0%", borderRadius: "10px", border: "2px solid white", bgcolor: colors.primary[800] }}></Box>
            <Box sx={{ width: "30px", height: "30px", position: "absolute", top: "50%", left: "0%", transform: "translate(-7px, -15px)", bgcolor: colors.primary[800] }}>
              <UploadIcon sx={{ width: "20px", height: "20px", position: "absolute", top: "50%", left: "50%", transform: "translate(-50%, -50%)" }} />
            </Box>
            <Box sx={{ width: { xs: "100%", sm: "50%", md: "50%", lg: "50%", xl: "50%" } }}>
              <Typography fontSize="1.2em" fontWeight="bold">Select & Upload</Typography>
              <Typography color={colors.neutral[500]}>Create a mapping project by choosing a reference atlas together with a pre-trained neural network (model) and upload your single-cell data.</Typography>
            </Box>
            <Box sx={{ width: { xs: "100%", sm: "50%", md: "50%", lg: "50%", xl: "50%" }, backgroundColor: "white", borderRadius: "20px" }}>
              <img style={{ width: "100%" }} src={graphic2} alt="Upload" />
            </Box>
          </Box>

          {/* PROCESSING */}
          <Box sx={{ position: "relative", gap: "1em", display: "flex", flexDirection: { xs: "column", sm: "row", md: "row", lg: "row", xl: "row" }, justifyContent: "space-between", alignItems: "center", width: "90%", margin: "auto", padding: "5% 10%" }}>
            <Box sx={{ position: "absolute", top: "0%", left: "7px", width: "2px", height: "100%", bgcolor: colors.neutral[700] }} />
            <Box sx={{ zIndex: "1", width: "30px", height: "30px", position: "absolute", top: "50%", left: "0%", transform: "translate(-7px, -15px)", bgcolor: colors.primary[800] }}>
              <AutorenewIcon sx={{ width: "20px", height: "20px", position: "absolute", top: "50%", left: "50%", transform: "translate(-50%, -50%)" }} />
            </Box>
            <Box sx={{ width: { xs: "100%", sm: "50%", md: "50%", lg: "50%", xl: "50%" } }}>
              <Typography fontSize="1.2em" fontWeight="bold">Processing</Typography>
              <Typography color={colors.neutral[500]}>ScArches will map your input data to the specified atlas using the selected machine learning model, while performing the cell-type classification.</Typography>
            </Box>
            <Box sx={{ width: { xs: "100%", sm: "50%", md: "50%", lg: "50%", xl: "50%" }, backgroundColor: "white", borderRadius: "20px" }}>
              <img style={{ width: "100%" }} src={graphic3} alt="Processing" />
            </Box>
          </Box>

          {/* CHECK RESULT */}
          <Box sx={{ position: "relative", gap: "1em", display: "flex", flexDirection: { xs: "column", sm: "row", md: "row", lg: "row", xl: "row" }, justifyContent: "space-between", alignItems: "center", width: "90%", margin: "auto", padding: "5% 10% 10% 10%" }}>
            <Box sx={{ position: "absolute", top: "0%", left: "7px", width: "2px", height: "calc(95% - 8px)", bgcolor: colors.neutral[700] }} />
            <Box sx={{ zIndex: "1", width: "16px", height: "16px", position: "absolute", top: "95%", left: "0%", transform: "translate(-0%, -100%)", bgcolor: colors.primary[800] }}>
              <CheckCircleOutlineIcon sx={{ width: "16px", height: "16px", color: colors.secondary1[400] }} />
            </Box>
            <Box sx={{ zIndex: "1", width: "30px", height: "30px", position: "absolute", top: "50%", left: "0%", transform: "translate(-7px, -15px)", bgcolor: colors.primary[800] }}>
              <AssignmentIcon sx={{ width: "20px", height: "20px", position: "absolute", top: "50%", left: "50%", transform: "translate(-50%, -50%)" }} />
            </Box>
            <Box sx={{ width: { xs: "100%", sm: "50%", md: "50%", lg: "50%", xl: "50%" } }}>
              <Typography fontSize="1.2em" fontWeight="bold">Check Results</Typography>
              <Typography color={colors.neutral[500]}>Your new mapping will be automaticly displayed for further analysis. </Typography>
            </Box>
            <Box sx={{ width: { xs: "100%", sm: "50%", md: "50%", lg: "50%", xl: "50%" }, backgroundColor: "white", borderRadius: "20px" }}>
              <img style={{ width: "100%" }} src={graphic4} alt="Check Results" />
            </Box>
          </Box>
        </Box>
        {/* CONTACT US */}
        <Box ref={contactUsBoxRef} sx={{ position: "relative", margin: "4em auto", position: "relative", width: "100%" }} >
          <Typography sx={{ textAlign: "center" }} fontSize="2em" fontWeight="bold">Contact Us</Typography>
          <Typography marginTop="1em" sx={{ textAlign: "center" }} fontSize="1em">Please message us in case you have any questions, feedback or collaboration-related inquiries concerning ArchMap.</Typography>
          {/* <Box sx={{
            width: "100%",
            margin: "2em auto 0em auto",
            padding: "2em 0em",
            boxShadow: "0px 4px 6px rgba(0, 0, 0, 0.10), 0px 0px 1px rgba(0, 0, 0, 0.20)",
            borderRadius: "10px"
          }}>
            <Stack sx={{ width: "80%", margin: "auto" }} direction="column" spacing={4}>
              <Stack direction="row" justifyContent="space-between" spacing={2}>
                <Input isRequired label="Email" />
                <Input isRequired label="First Name" />
                <Input isRequired label="Last Name" />
              </Stack>
              <Input isRequired label="Subject" maxLength={50}/>
              <Input isRequired label="Message" multiline maxLength={1000} />
              <Box sx={{ display: "flex", flexDirection: "row", justifyContent: "center" }}>
                <CustomButton>Send</CustomButton>
              </Box>
            </Stack>
          </Box> */}
          <ContactForm />
        </Box>
      </Box>

      {/* FOOTER */}
      <Footer />
    </Box>
  )
}

export default Home