import CircleIcon from '@mui/icons-material/Circle';
import { Box, CardActionArea, Stack } from "@mui/material";
import Typography from '@mui/material/Typography';
import { useEffect, useState } from 'react';
import { colors } from 'shared/theme/colors';

/** 
 * TabCard for TabGroup Component
 * @param width
 * @param height
 * @status status of file upload
 * @fileName filename
 * @handleOnClick executed function on card click
*/
export const TabCard = ({ width, height, status, fileName, handleOnClick }) => {
    const [color, setColor] = useState(colors.success.main)
    // TODO upload percentage?

    // useEffect(() => {
        // TODO this part is not working
        // status === 'DONE' ? 
        // setColor(colors.success.main) :
        // status === 'IN PROGRESS' ?
        // setColor(colors.warning.main) :
        // setColor(colors.red.main)
    // }, [color])

    return (
        <CardActionArea disableTouchRipple onClick={handleOnClick}>
            <Box sx={{
                    width: {width},
                    height: {height},
                    backgroundColor: 'white',
                    borderRadius: "0.625rem",
                    marginTop: '0.5em',
                    ":hover": {
                        color: colors.primary,
                        ':hover': { backgroundColor: colors.primary[300], transition: '0.4s', color: 'white' },
                        ':focus': { backgroundColor: colors.primary[300], transition: '0.4s', color: 'white' },
                        ':disabled': { backgroundColor: '#EBEFFF', transition: '0.4s', color: colors.primary[600] }
                    }
                }}
            >
                <Box sx={{ boxShadow: "0px 0px 2px rgba(0,0,0, 0.15)", p: "0.5em", borderRadius: "0.625rem" }}>
                    <Stack direction="row" spacing={4}>
                        <CircleIcon sx={{ fontSize: 18, marginLeft: '3%', marginTop: '0.2em', color }} />
                        <Typography variant='body1'>{fileName}</Typography>
                    </Stack>
                </Box>
            </Box>
        </CardActionArea>
    );
}