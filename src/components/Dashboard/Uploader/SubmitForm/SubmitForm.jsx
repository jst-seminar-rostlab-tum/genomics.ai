import styles from "./submitform.module.css";
import { useState } from "react";
import { Button } from "@mui/material";
import UploadField from "../UploadField/UploadField";
import ProgressBar from "../ProgressBar/ProgressBar";

function SubmitForm(){
    const [fileAvailable, setFileAvailable] = useState(false);
    const [selectedFile, setSelectedFile] = useState();
    //verification complete might not be a necessary constant
    const [verificationComplete, setVerificationComplete] = useState(true);
    const [fileBeingUploaded, setFileBeingUploaded] = useState(false);
    const [fileBeingVerified, setFileBeingVerified] = useState(false);
    const [uploadFieldDisabled, setUploadFieldDisabled] = useState(false);
    const backendUrl = 'https://custom-helix-329116.ey.r.appspot.com/hello';

    let button;
    //type of submit button, calls the submission handler when available
    if(fileAvailable && !fileBeingUploaded){
        button = <Button 
        onClick={()=> {
            setFileBeingUploaded(true);
            handleSubmission();
        }}>Submit</Button>;
    }
    else{
        button = <Button disabled>Submit</Button>;
    }

    //button status during file verification
    if(fileBeingVerified){
        button = <Button disabled>Verifying</Button>
    }
    
    //handling of the selected file for upload
    const changeHandler = (event) => {
        //this function returns true if the file is valid for submission
        function inputValidation(){
            let validFile = true;
            //todo: call api for input validation
            //block the sequence until the verification is done

            return validFile;
        }

        setFileBeingVerified(true);
        setUploadFieldDisabled(true);
        setSelectedFile(event.target.files[0]);
        //setFileBeingVerified(false); 

        if(inputValidation()){
            setFileAvailable(true);
            setFileBeingVerified(false);
        }else{
            setFileAvailable(false);
            setFileBeingVerified(false);
            setUploadFieldDisabled(false);
            //todo, return error message
        }
    }

    //uploads the file upon submission of the file
    const handleSubmission = () => {
        //example fetch
        fetch(backendUrl, {
            method: 'POST',
            body: selectedFile
        })
        //call the multipart upload here
    }

    return (
        <div className={styles.SubmitForm}>
            <UploadField 
            disabled={uploadFieldDisabled} 
            selectedFile={selectedFile}
            onChange={changeHandler}
            /> 
            <div className={styles.progressTitleUpdate}>
                { fileBeingUploaded ? "uploading..." : null}
            </div>
            <div className={styles.flexContainer}>
                {fileBeingUploaded ? <ProgressBar value={submissionProgress}/> : null}
            </div>
            {/* <div>
                <Snackbar open={open} autoHideDuration={6000} onClose={handleClose}>
                    <Alert onClose={handleClose} severity="success" sx={{ width: '100%' }}>
                    This is a success message!
                    </Alert>
                </Snackbar>
            </div> */}
            <div>
                {button}
            </div>
        </div>
    )

    //upload field
    //submit button
}



//not yet available in the progress bar as a prop
function submissionProgress(){
    
}


export default SubmitForm;