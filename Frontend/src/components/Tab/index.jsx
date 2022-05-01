import React, { useRef, useState, useEffect } from 'react';
import { useLocation, Link } from 'react-router-dom';

import { Tab, Tabs, Box } from '@mui/material';

/**
 * Styled Tab
 * packaged the Tab from the MUI and add style to fit our mock
 *
 * need a boolean parameter
 * - darkBackground
 * in order to render the Tab accordingly
 *
 * and a label, the content
 *
 * see the example StyledTab in TabGroup below for how to combine the Tab with a Link from react-router-dom
 */
export function StyledTab(props) {
  const {
    darkBackground, label, component, to, ...rest
  } = props;

  return (
    <Tab
      {...rest}
      disableRipple
      label={label}
      component={component}
      to={to}
      sx={{
        textTransform: 'none',
        fontWeight: 'bold',
        color: darkBackground ? 'white' : 'black',
        opacity: 1,
        '&.Mui-selected': {
          fontWeight: 'bold',
          color: darkBackground ? 'white' : 'black',
          opacity: 1,
        },
      }}
    />
  );
}

/**
 * Tab Group, the collection of Tabs
 *
 * !!! important
 * the value state is now designed to be organizd by it's parent to make it more flexble
 *
 * like the StyledTab, it also need a boolean parameter
 * - darkBackground
 * to render the Tab accordingly
 *
 * it also need an array of objects containing labels and paths / additionalContent
 *
 * for example:
 * [
 *     {
 *         label: "tab 1",
 *         path: "/tab_1"
 *     },
 *     {
 *         label: "tab 2",
 *         path: "/tab_2"
 *     },
 *     {
 *         label: "tab 3",
 *         path: "/tab_3"
 *     }
 * ]
 *
 * or:
 * [
 *     {
 *         label: "tab 1",
 *         additionalContent: (<Box>something</Box>)
 *     },
 *     {
 *         label: "tab 2",
 *         additionalContent: (<Box>something else</Box>)
 *     }
 * ]
 *
 * the default width and height is 100% of the containing block
 * but it is also possible to customize it
 *
 */
export function TabGroup(props) {
  const position = useLocation();

  const {
    value, setValue, darkBackground, tabsInfo, width = '100%', height = '100%',
  } = props;

  const tabsRef = useRef();
  const [tabsHeight, setTabsHeight] = useState(0);

  useEffect(() => {
    setTabsHeight(tabsRef.current.clientHeight);
  }, []);

  return (
    <Box sx={{ width, height, position: 'relative' }}>
      <Box ref={tabsRef}>
        <Tabs
          value={value}
          onChange={(_, newValue) => setValue(newValue)}
          sx={{
            '& .MuiTabs-indicator': {
              backgroundColor: darkBackground ? 'white' : 'black',
            },
          }}
        >
          {
            tabsInfo.map((tabInfo, index) => (
              <StyledTab
                key={tabInfo.label}
                value={index}
                label={tabInfo.label}
                component={Link}
                to={tabInfo.path ? tabInfo.path : position.pathname}
                darkBackground={darkBackground}
              />
            ))
                    }
        </Tabs>
      </Box>

      {
                tabsInfo[value].additionalContent
                && (
                <Box sx={{ width: '100%', height: `calc(100% - ${tabsHeight}px)`, overflowY: 'scroll' }}>
                  {tabsInfo[value].additionalContent}
                </Box>
                )
            }

    </Box>
  );
}
