import {Box, Button} from "@mui/material";
import styles from "./topbar.module.css"
import {Link} from "react-router-dom";


function TopBar() {
    return (
        <div className={styles.container}>
            <Box mx={1}>
                <Button variant="contained" component={Link} to={"/"}>
                    <span className="material-icons">home</span>&nbsp;
                    Page One
                </Button>
            </Box>
            <Box mx={1}>
                <Button variant="contained" component={Link} to={"/page-two"}>
                    Page Two
                </Button>
            </Box>
        </div>
    );
}

export default TopBar;