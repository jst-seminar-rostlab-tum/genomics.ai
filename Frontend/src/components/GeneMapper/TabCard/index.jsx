import CircleIcon from '@mui/icons-material/Circle';
import { Box, CardActionArea, Stack } from "@mui/material";
import Typography from '@mui/material/Typography';
import { useEffect, useState } from 'react';
import { colors } from 'shared/theme/colors';

/** 
 * TabCard for TabGroup Component / List of cards with select feature used in FileUpload page
 * @param width
 * @param height
 * @data object containing information to be displayed on the card
 * @handleOnClick executed function on card click
 * @selected boolean parameter to indicate if a card element is selected
*/
export const TabCard = ({ width, height, data, handleOnClick, selected}) => {
    // commented out for (possibly) later use
    // const [color, setColor] = useState(colors.success.main);

    // useEffect(() => {
        // TODO this part is not working
        // status === 'DONE' ? 
        // setColor(colors.success.main) :
        // status === 'IN PROGRESS' ?
        // setColor(colors.warning.main) :
        // setColor(colors.red.main)
    // }, [color])

    return (
        <Box onClick={handleOnClick} sx={{ width: {width}, height: {height}, backgroundColor: 'white', borderRadius: "0.625rem", marginTop: '0.5em' }}>
            <Box 
                sx={ selected ? {
                    boxShadow: "0px 0px 2px rgba(0,0,0, 0.15)", 
                    p: "0.5em", 
                    borderRadius: "0.625rem",
                    backgroundColor: colors.primary[300],
                    color: 'white'
                }
                : { 
                boxShadow: "0px 0px 2px rgba(0,0,0, 0.15)", 
                p: "0.5em", 
                borderRadius: "0.625rem",
                backgroundColor: 'white',
                ":hover": {
                    color: colors.primary,
                    ':hover': { backgroundColor: colors.primary[300], transition: '0.4s', color: 'white' },
                    ':focus': { backgroundColor: colors.primary[300], transition: '0.4s', color: 'white' },
                    ':disabled': { backgroundColor: '#EBEFFF', transition: '0.4s', color: colors.primary[600] }
                }}}
            >
                <Stack p='0.1em' pl='0.3em'>
                    {/* <CircleIcon sx={{ fontSize: 18, marginLeft: '3%', marginTop: '0.2em', color }} /> */}
                    <Typography variant='body1'>
                        {data.name}
                    </Typography>
                    <Typography variant='caption'>
                        {data.visibility}
                    </Typography>
                </Stack>
            </Box>
        </Box>
    );
}